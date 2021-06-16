
#ifdef MACOSX
// This is to overcome a possible bug in Apple's GCC.
#define isnan(x) (x==FP_NAN)
#endif

/**Just scan match every single particle.
If the scan matching fails, the particle gets a default likelihood.*/
/// 为每一个粒子都计算扫描匹配，并计算每个粒子的权重
inline void GridSlamProcessor::scanMatch(const double* plainReading)
{
  // sample a new pose from each scan in the reference
  /*每个粒子都要进行scan-match*/
  double sumScore=0;
  
  /*遍历每一个粒子，为每一个粒子都计算扫描匹配，并计算权重*/
  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
  {
    /*对当前粒子进行扫描匹配，计算最优位姿*/
    OrientedPoint corrected;
    double score, l, s;
    score=m_matcher.optimize(corrected, it->map, it->pose, plainReading);
    
    /*矫正成功则更新位姿*/
    //    it->pose=corrected;
    if (score>m_minimumScore)
    {
      it->pose=corrected;
    } 
    else 
    {
        if (m_infoStream)
        {
          m_infoStream << "Scan Matching Failed, using odometry. Likelihood=" << l <<std::endl;
          m_infoStream << "lp:" << m_lastPartPose.x << " "  << m_lastPartPose.y << " "<< m_lastPartPose.theta <<std::endl;
          m_infoStream << "op:" << m_odoPose.x << " " << m_odoPose.y << " "<< m_odoPose.theta <<std::endl;
        }
    }

    // 粒子的最优位姿计算了之后，重新计算粒子的权重
    // 相当于粒子滤波器中的观测步骤，计算p(z|x,m)，粒子的权重由粒子的似然来表示
    m_matcher.likelihoodAndScore(s, l, it->map, it->pose, plainReading);
    sumScore+=score;
    it->weight+=l;
    it->weightSum+=l;

    //set up the selective copy of the active area
    //by detaching the areas that will be updated
    m_matcher.invalidateActiveArea();
    m_matcher.computeActiveArea(it->map, it->pose, plainReading);
  }
  if (m_infoStream)
    m_infoStream << "Average Scan Matching Score=" << sumScore/m_particles.size() << std::endl;	
}

/// 归一化粒子权重，并计算权重相似度Neff（用于判断是否需要重采样）
/// Neff​越大说明粒子越收敛，越小说明越分散
inline void GridSlamProcessor::normalize()
{
  //normalize the log m_weights
  double gain=1./(m_obsSigmaGain*m_particles.size());
  
  /*求所有粒子中的最大的权重*/
  double lmax= -std::numeric_limits<double>::max();
  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
  {
    lmax=it->weight>lmax?it->weight:lmax;
  }
  //cout << "!!!!!!!!!!! maxwaight= "<< lmax << endl;
  
  /*以最大权重为中心的高斯分布*/
  m_weights.clear();
  double wcum=0;
  m_neff=0;
  for (std::vector<Particle>::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
  {
    m_weights.push_back(exp(gain*(it->weight-lmax)));
    wcum+=m_weights.back();
    //cout << "l=" << it->weight<< endl;
  }
  
  /*归一化权重*/
  m_neff=0; // 权重相似度
  for (std::vector<double>::iterator it=m_weights.begin(); it!=m_weights.end(); it++)
  {
    *it=*it/wcum;

    double w=*it;
    m_neff+=w*w;
  }
  m_neff=1./m_neff;
  
}

/// 粒子滤波器重采样
inline bool GridSlamProcessor::resample(const double* plainReading, int adaptSize, const RangeReading* reading)
{
  
  bool hasResampled = false;
  
  /*备份老的粒子的轨迹  即保留叶子节点 在增加新节点的时候使用*/
  TNodeVector oldGeneration;
  for (unsigned int i=0; i<m_particles.size(); i++)
  {
    oldGeneration.push_back(m_particles[i].node);
  }
  
  if (m_neff<m_resampleThreshold*m_particles.size())
  {	
    /// -------------------------需要重采样------------------------
    
    if (m_infoStream)
      m_infoStream  << "*************RESAMPLE***************" << std::endl;
    
    // 采取重采样方法决定 哪些粒子会保留
    // 保留的粒子会返回下标 里面的下标可能会重复 因为有些粒子会重复采样，而另外的一些粒子会消失掉
    uniform_resampler<double, double> resampler;
    m_indexes=resampler.resampleIndexes(m_weights, adaptSize);
    
    if (m_outputStream.is_open())
    {
      m_outputStream << "RESAMPLE "<< m_indexes.size() << " ";
      for (std::vector<unsigned int>::const_iterator it=m_indexes.begin(); it!=m_indexes.end(); it++)
      {
          m_outputStream << *it <<  " ";
      }
      m_outputStream << std::endl;
    }
    
    onResampleUpdate();
    //BEGIN: BUILDING TREE
    
    // 重采样之后的粒子
    ParticleVector temp;
    unsigned int j=0;
    
    // 要删除的粒子下标
    std::vector<unsigned int> deletedParticles; //this is for deleteing the particles which have been resampled away.
    
    //		cerr << "Existing Nodes:" ;
    // 枚举每一个要被保留的粒子
    for (unsigned int i=0; i<m_indexes.size(); i++)
    {
      //			cerr << " " << m_indexes[i];
      // 统计要被删除的粒子
      while(j<m_indexes[i])
      {
        deletedParticles.push_back(j);
        j++;
      }
      if (j==m_indexes[i])
          j++;
      
      // 得到当前的保留的粒子
      Particle & p=m_particles[m_indexes[i]];
      
      // 每一个需要保留下来的粒子都需要在路径中增加一个新的节点
      TNode* node=0;
      TNode* oldNode=oldGeneration[m_indexes[i]];
      //			cerr << i << "->" << m_indexes[i] << "B("<<oldNode->childs <<") ";
      
      // 创建一个新的节点 改节点的父节点为oldNode
      node=new	TNode(p.pose, 0, oldNode, 0);
      //node->reading=0;
      node->reading=reading;
      //			cerr << "A("<<node->parent->childs <<") " <<endl;
      
      // 这个要保留下来的粒子，要保留的粒子的下标为m_indexs
      temp.push_back(p);
      temp.back().node=node;
      temp.back().previousIndex=m_indexes[i];
    }
    
    while(j<m_indexes.size())
    {
      deletedParticles.push_back(j);
      j++;
    }
    
    //		cerr << endl;
    std::cerr <<  "Deleting Nodes:";
    // 把要删除的粒子的Node都删除掉，Node表示轨迹的起点(最新的点)
    for (unsigned int i=0; i<deletedParticles.size(); i++)
    {
      std::cerr <<" " << deletedParticles[i];
      delete m_particles[deletedParticles[i]].node;
      m_particles[deletedParticles[i]].node=0;
    }
    std::cerr  << " Done" <<std::endl;
    
    //END: BUILDING TREE
    std::cerr << "Deleting old particles..." ;
    // 清除全部的粒子 然后从tmp中读取保留下来的粒子，并更新粒子的地图
    m_particles.clear();
    std::cerr << "Done" << std::endl;
    std::cerr << "Copying Particles and  Registering  scans...";
    for (ParticleVector::iterator it=temp.begin(); it!=temp.end(); it++)
    {
      it->setWeight(0);
      m_matcher.invalidateActiveArea();
      m_matcher.registerScan(it->map, it->pose, plainReading);
      m_particles.push_back(*it);
    }
    std::cerr  << " Done" <<std::endl;
    hasResampled = true;
  } 
  else 
  {
    /// -----------------------不需要重采样-------------------------
      
    /// 保留所有的粒子，为每个粒子创建新的节点，并更新所有粒子的地图
    int index=0;
    std::cerr << "Registering Scans:";
    TNodeVector::iterator node_it=oldGeneration.begin();
    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
    {
      //create a new node in the particle tree and add it to the old tree
      //BEGIN: BUILDING TREE  
      // 为当前粒子创建 新的节点
      TNode* node=0;
      node=new TNode(it->pose, 0.0, *node_it, 0);
      
      //node->reading=0;
      node->reading=reading;
      it->node=node;

      //END: BUILDING TREE
      // 更新粒子的地图
      m_matcher.invalidateActiveArea();
      m_matcher.registerScan(it->map, it->pose, plainReading);
      it->previousIndex=index;
      index++;
      node_it++;
      
    }
    std::cerr  << "Done" <<std::endl;
    
  }
  //END: BUILDING TREE
  
  return hasResampled;
}
