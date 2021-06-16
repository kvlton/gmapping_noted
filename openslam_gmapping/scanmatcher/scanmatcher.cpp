#include <cstring>
#include <limits>
#include <list>
#include <iostream>

#include "gmapping/scanmatcher/scanmatcher.h"
#include "gmapping/scanmatcher/gridlinetraversal.h"
//#define GENERATE_MAPS

namespace GMapping {

using namespace std;

const double ScanMatcher::nullLikelihood=-.5;

ScanMatcher::ScanMatcher(): m_laserPose(0,0,0){
	//m_laserAngles=0;
	m_laserBeams=0;
	m_optRecursiveIterations=3;
	m_activeAreaComputed=false;

	// This  are the dafault settings for a grid map of 5 cm
	m_llsamplerange=0.01;
	m_llsamplestep=0.01;
	m_lasamplerange=0.005;
	m_lasamplestep=0.005;
	m_enlargeStep=10.;
	m_fullnessThreshold=0.1;
	m_angularOdometryReliability=0.;
	m_linearOdometryReliability=0.;
	m_freeCellRatio=sqrt(2.);
	m_initialBeamsSkip=0;
	
/*	
	// This  are the dafault settings for a grid map of 10 cm
	m_llsamplerange=0.1;
	m_llsamplestep=0.1;
	m_lasamplerange=0.02;
	m_lasamplestep=0.01;
*/	
	// This  are the dafault settings for a grid map of 20/25 cm
/*
	m_llsamplerange=0.2;
	m_llsamplestep=0.1;
	m_lasamplerange=0.02;
	m_lasamplestep=0.01;
	m_generateMap=false;
*/

   m_linePoints = new IntPoint[20000];
}

ScanMatcher::~ScanMatcher(){
	delete [] m_linePoints;
}

void ScanMatcher::invalidateActiveArea(){
	m_activeAreaComputed=false;
}

/*
void ScanMatcher::computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	if (m_activeAreaComputed)
		return;
	HierarchicalArray2D<PointAccumulator>::PointSet activeArea;
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
	const double * angle=m_laserAngles;
	for (const double* r=readings; r<readings+m_laserBeams; r++, angle++)
		if (m_generateMap){
			double d=*r;
			if (d>m_laserMaxRange)
				continue;
			if (d>m_usableRange)
				d=m_usableRange;
			
			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p1=map.world2map(phit);
			
			d+=map.getDelta();
			//Point phit2=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			//IntPoint p2=map.world2map(phit2);
			IntPoint linePoints[20000] ;
			GridLineTraversalLine line;
			line.points=linePoints;
			//GridLineTraversal::gridLine(p0, p2, &line);
			GridLineTraversal::gridLine(p0, p1, &line);
			for (int i=0; i<line.num_points-1; i++){
				activeArea.insert(map.storage().patchIndexes(linePoints[i]));
			}
			if (d<=m_usableRange){
				activeArea.insert(map.storage().patchIndexes(p1));
				//activeArea.insert(map.storage().patchIndexes(p2));
			}
		} else {
			if (*r>m_laserMaxRange||*r>m_usableRange) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			IntPoint cp=map.storage().patchIndexes(p1);
			assert(cp.x>=0 && cp.y>=0);
			activeArea.insert(cp);
			
		}
	//this allocates the unallocated cells in the active area of the map
	//cout << "activeArea::size() " << activeArea.size() << endl;
	map.storage().setActiveArea(activeArea, true);
	m_activeAreaComputed=true;
}
*/


/// 计算有效区域，参数：概率栅格图 机器人位置 激光雷达数据
void ScanMatcher::computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings)
{
    /*已经计算过了，则没必要计算了*/
    if (m_activeAreaComputed)
		return;
    
    /*把激光雷达的位置转换到地图坐标系*/
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
	
    ///------------------------------扩充地图范围-------------------------------
    
    /*地图的范围*/
	Point min(map.map2world(0,0));
	Point max(map.map2world(map.getMapSizeX()-1,map.getMapSizeY()-1));
	       
    /*扩展地图的范围*/
	if (lp.x<min.x) min.x=lp.x;
	if (lp.y<min.y) min.y=lp.y;
	if (lp.x>max.x) max.x=lp.x;
	if (lp.y>max.y) max.y=lp.y;
	
	/*determine the size of the area*/
    /*遍历激光数据，扩展地图的范围*/
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
    {
        /*去除不合理的值*/
		if (*r>m_laserMaxRange||*r==0.0||isnan(*r)) continue;
        
        /*截断*/
		double d=*r>m_usableRange?m_usableRange:*r;
        
        /*被激光击中的位置*/
		Point phit=lp;
		phit.x+=d*cos(lp.theta+*angle);
		phit.y+=d*sin(lp.theta+*angle);
        
        /*扩充地图范围*/
		if (phit.x<min.x) min.x=phit.x;
		if (phit.y<min.y) min.y=phit.y;
		if (phit.x>max.x) max.x=phit.x;
		if (phit.y>max.y) max.y=phit.y;
	}
	//min=min-Point(map.getDelta(),map.getDelta());
	//max=max+Point(map.getDelta(),map.getDelta());
	
    /*如果地图需要扩展*/
	if ( !map.isInside(min)	|| !map.isInside(max))
    {
        /*得到目前地图的大小*/
		Point lmin(map.map2world(0,0));
		Point lmax(map.map2world(map.getMapSizeX()-1,map.getMapSizeY()-1));
		//cerr << "CURRENT MAP " << lmin.x << " " << lmin.y << " " << lmax.x << " " << lmax.y << endl;
		//cerr << "BOUNDARY OVERRIDE " << min.x << " " << min.y << " " << max.x << " " << max.y << endl;
        
        /*地图扩充*/
		min.x=( min.x >= lmin.x )? lmin.x: min.x-m_enlargeStep;
		max.x=( max.x <= lmax.x )? lmax.x: max.x+m_enlargeStep;
		min.y=( min.y >= lmin.y )? lmin.y: min.y-m_enlargeStep;
		max.y=( max.y <= lmax.y )? lmax.y: max.y+m_enlargeStep;
		map.resize(min.x, min.y, max.x, max.y);
		//cerr << "RESIZE " << min.x << " " << min.y << " " << max.x << " " << max.y << endl;
	}
    
    ///-------------------------------分配有效区域------------------------------------
	
	HierarchicalArray2D<PointAccumulator>::PointSet activeArea;
	/*allocate the active area*/
    /*遍历激光数据，分配有效区域*/
	angle=m_laserAngles+m_initialBeamsSkip;
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
		if (m_generateMap)
        {
            /// 需要生成地图
            
            // 排除错误的激光点
			double d=*r;
			if (d>m_laserMaxRange||d==0.0||isnan(d))
				continue;
			if (d>m_usableRange)
				d=m_usableRange;
            
            /*激光束的起点和终点*/
			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p0=map.world2map(lp);   // 起点
			IntPoint p1=map.world2map(phit); // 终点
			
			//IntPoint linePoints[20000] ;
            /*bresenham算法来计算激光起点到终点要经过的路径*/
			GridLineTraversalLine line;
			line.points=m_linePoints;
			GridLineTraversal::gridLine(p0, p1, &line);
            
            /*更新地图 把画线算法计算出来的值都算进去*/
			for (int i=0; i<line.num_points-1; i++)
            {
				assert(map.isInside(m_linePoints[i]));
				activeArea.insert(map.storage().patchIndexes(m_linePoints[i]));
				assert(m_linePoints[i].x>=0 && m_linePoints[i].y>=0);
			}
            
            /*如果终点在可用范围，也加入有效范围*/
			if (d<m_usableRange)
            {
				IntPoint cp=map.storage().patchIndexes(p1);
				assert(cp.x>=0 && cp.y>=0);
				activeArea.insert(cp);
			}
		} 
        else 
        {
            /// 不需要生成地图
            
			if (*r>m_laserMaxRange||*r>m_usableRange||*r==0.0||isnan(*r)) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			IntPoint cp=map.storage().patchIndexes(p1);
			assert(cp.x>=0 && cp.y>=0);
			activeArea.insert(cp);
		}
	
	//this allocates the unallocated cells in the active area of the map
	//cout << "activeArea::size() " << activeArea.size() << endl;
/*	
	cerr << "ActiveArea=";
	for (HierarchicalArray2D<PointAccumulator>::PointSet::const_iterator it=activeArea.begin(); it!= activeArea.end(); it++){
		cerr << "(" << it->x <<"," << it->y << ") ";
	}
	cerr << endl;
*/		
	map.storage().setActiveArea(activeArea, true);
	m_activeAreaComputed=true;
}

/// 激光数据注册，参数：概率栅格图 机器人位置 激光雷达数据
double ScanMatcher::registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings)
{
	if (!m_activeAreaComputed)
		computeActiveArea(map, p, readings);
		
	//this operation replicates the cells that will be changed in the registration operation
	map.storage().allocActiveArea();
	
    /*把激光雷达的位置转换到地图坐标系*/
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp); // 激光雷达在地图坐标系的位置，即所有激光束的起点
	
	/*遍历激光数据，更新栅格地图*/
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	double esum=0; // 表示总体的熵的增加（或者减少）
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
		if (m_generateMap)
        {
            /// 需要生成地图
            
            // 排除错误的激光点
			double d=*r;
			if (d>m_laserMaxRange||d==0.0||isnan(d))
				continue;
			if (d>m_usableRange)
				d=m_usableRange;
            
            /*bresenham算法来计算激光起点到终点要经过的路径*/
			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p1=map.world2map(phit);
			//IntPoint linePoints[20000] ;
			GridLineTraversalLine line;
			line.points=m_linePoints;
			GridLineTraversal::gridLine(p0, p1, &line);
            
            /*更新空闲位置*/
			for (int i=0; i<line.num_points-1; i++)
            {
				PointAccumulator& cell=map.cell(line.points[i]);
				double e=-cell.entropy();
				cell.update(false, Point(0,0)); // free
				e+=cell.entropy();
				esum+=e;
			}
            
            /*更新被击中的位置，用来标记障碍物*/
			if (d<m_usableRange)
            {
				double e=-map.cell(p1).entropy();
				map.cell(p1).update(true, phit); // occupy
				e+=map.cell(p1).entropy();
				esum+=e;
			}
		} 
        else 
        {
            /// 不需要生成地图
            
			if (*r>m_laserMaxRange||*r>m_usableRange||*r==0.0||isnan(*r)) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			map.cell(p1).update(true,phit);
		}
	//cout  << "informationGain=" << -esum << endl;
	return esum;
}

/*
void ScanMatcher::registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	if (!m_activeAreaComputed)
		computeActiveArea(map, p, readings);
		
	//this operation replicates the cells that will be changed in the registration operation
	map.storage().allocActiveArea();
	
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
	const double * angle=m_laserAngles;
	for (const double* r=readings; r<readings+m_laserBeams; r++, angle++)
		if (m_generateMap){	
			double d=*r;
			if (d>m_laserMaxRange)
				continue;
			if (d>m_usableRange)
				d=m_usableRange;
			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p1=map.world2map(phit);
			
			IntPoint linePoints[20000] ;
			GridLineTraversalLine line;
			line.points=linePoints;
			GridLineTraversal::gridLine(p0, p1, &line);
			for (int i=0; i<line.num_points-1; i++){
				IntPoint ci=map.storage().patchIndexes(line.points[i]);
				if (map.storage().getActiveArea().find(ci)==map.storage().getActiveArea().end())
					cerr << "BIG ERROR" <<endl;
				map.cell(line.points[i]).update(false, Point(0,0));
			}
			if (d<=m_usableRange){
				
				map.cell(p1).update(true,phit);
			}
		} else {
			if (*r>m_laserMaxRange||*r>m_usableRange) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			map.cell(phit).update(true,phit);
		}
}

*/

double ScanMatcher::icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	double currentScore;
	double sc=score(map, init, readings);;
	OrientedPoint start=init;
	pnew=init;
	int iterations=0;
	do{
		currentScore=sc;
		sc=icpStep(pnew, map, start, readings);
		//cerr << "pstart=" << start.x << " " <<start.y << " " << start.theta << endl;
		//cerr << "pret=" << pnew.x << " " <<pnew.y << " " << pnew.theta << endl;
		start=pnew;
		iterations++;
	} while (sc>currentScore);
	cerr << "i="<< iterations << endl;
	return currentScore;
}

/// 根据地图、激光数据、位姿迭代求解一个最优的新的位姿出来
double ScanMatcher::optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const
{
	double bestScore=-1;
    
    /*计算当前位置的得分*/
	OrientedPoint currentPose=init;
	double currentScore=score(map, currentPose, readings);
    
    /*搜索的步长（距离步长、角度步长），会逐渐缩短*/
	double adelta=m_optAngularDelta, ldelta=m_optLinearDelta;
    
	unsigned int refinement=0; // 大迭代的次数，每缩短一次步长表示进行了一次大迭代
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done}; // 用于小迭代的枚举类型
/*	cout << __PRETTY_FUNCTION__<<  " readings: ";
	for (int i=0; i<m_laserBeams; i++){
		cout << readings[i] << " ";
	}
	cout << endl;
*/	int c_iterations=0; // 小迭代的次数
	do{
        /*如果这一次迭代(currentScore)算出来的比上一次(bestScore)差，则有可能是走太多了，要减少搜索步长 这个策略跟LM有点像*/
		if (bestScore>=currentScore)
        {
            /*缩短步长，进行下一次大迭代*/
			refinement++;
			adelta*=.5;
			ldelta*=.5;
		}
        
		bestScore=currentScore;
//		cout <<"score="<< currentScore << " refinement=" << refinement;
//		cout <<  "pose=" << currentPose.x  << " " << currentPose.y << " " << currentPose.theta << endl;
		OrientedPoint bestLocalPose=currentPose;
		OrientedPoint localPose=currentPose;

        /*把8个方向都搜索一次  得到这8个方向里面最好的一个位姿和对应的得分*/
		Move move=Front;
		do {
			localPose=currentPose;
			switch(move){
				case Front:
					localPose.x+=ldelta;
					move=Back;
					break;
				case Back:
					localPose.x-=ldelta;
					move=Left;
					break;
				case Left:
					localPose.y-=ldelta;
					move=Right;
					break;
				case Right:
					localPose.y+=ldelta;
					move=TurnLeft;
					break;
				case TurnLeft:
					localPose.theta+=adelta;
					move=TurnRight;
					break;
				case TurnRight:
					localPose.theta-=adelta;
					move=Done;
					break;
				default:;
			}
			
            //计算当前的位姿和初始位姿的区别 区别越大增益越小
			double odo_gain=1;
			if (m_angularOdometryReliability>0.){
				double dth=init.theta-localPose.theta; 	dth=atan2(sin(dth), cos(dth)); 	dth*=dth;
				odo_gain*=exp(-m_angularOdometryReliability*dth);
			}
			if (m_linearOdometryReliability>0.){
				double dx=init.x-localPose.x;
				double dy=init.y-localPose.y;
				double drho=dx*dx+dy*dy;
				odo_gain*=exp(-m_linearOdometryReliability*drho);
			}
            
            /*计算得分=增益*score*/
			double localScore=odo_gain*score(map, localPose, readings);
			
            /*如果得分更好，则更新*/
			if (localScore>currentScore){
				currentScore=localScore;
				bestLocalPose=localPose;
			}
			c_iterations++; // 进行了一次小迭代
		} while(move!=Done);
        
		currentPose=bestLocalPose;
//		cout << "currentScore=" << currentScore<< endl;
		//here we look for the best move;
	}while (currentScore>bestScore || refinement<m_optRecursiveIterations); // 直到得分不再增加，并且已经进行了3次大迭代
	//cout << __PRETTY_FUNCTION__ << "bestScore=" << bestScore<< endl;
	//cout << __PRETTY_FUNCTION__ << "iterations=" << c_iterations<< endl;
	
    /*返回最优位置和得分*/
    pnew=currentPose;
	return bestScore;
}

struct ScoredMove{
	OrientedPoint pose;
	double score;
	double likelihood;
};

typedef std::list<ScoredMove> ScoredMoveList;

double ScanMatcher::optimize(OrientedPoint& _mean, ScanMatcher::CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	ScoredMoveList moveList;
	double bestScore=-1;
	OrientedPoint currentPose=init;
	ScoredMove sm={currentPose,0,0};
	unsigned int matched=likelihoodAndScore(sm.score, sm.likelihood, map, currentPose, readings);
	double currentScore=sm.score;
	moveList.push_back(sm);
	double adelta=m_optAngularDelta, ldelta=m_optLinearDelta;
	unsigned int refinement=0;
	int count=0;
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done};
	do{
		if (bestScore>=currentScore){
			refinement++;
			adelta*=.5;
			ldelta*=.5;
		}
		bestScore=currentScore;
//		cout <<"score="<< currentScore << " refinement=" << refinement;
//		cout <<  "pose=" << currentPose.x  << " " << currentPose.y << " " << currentPose.theta << endl;
		OrientedPoint bestLocalPose=currentPose;
		OrientedPoint localPose=currentPose;

		Move move=Front;
		do {
			localPose=currentPose;
			switch(move){
				case Front:
					localPose.x+=ldelta;
					move=Back;
					break;
				case Back:
					localPose.x-=ldelta;
					move=Left;
					break;
				case Left:
					localPose.y-=ldelta;
					move=Right;
					break;
				case Right:
					localPose.y+=ldelta;
					move=TurnLeft;
					break;
				case TurnLeft:
					localPose.theta+=adelta;
					move=TurnRight;
					break;
				case TurnRight:
					localPose.theta-=adelta;
					move=Done;
					break;
				default:;
			}
			double localScore, localLikelihood;
			
			double odo_gain=1;
			if (m_angularOdometryReliability>0.){
				double dth=init.theta-localPose.theta; 	dth=atan2(sin(dth), cos(dth)); 	dth*=dth;
				odo_gain*=exp(-m_angularOdometryReliability*dth);
			}
			if (m_linearOdometryReliability>0.){
				double dx=init.x-localPose.x;
				double dy=init.y-localPose.y;
				double drho=dx*dx+dy*dy;
				odo_gain*=exp(-m_linearOdometryReliability*drho);
			}
			localScore=odo_gain*score(map, localPose, readings);
			//update the score
			count++;
			matched=likelihoodAndScore(localScore, localLikelihood, map, localPose, readings);
			if (localScore>currentScore){
				currentScore=localScore;
				bestLocalPose=localPose;
			}
			sm.score=localScore;
			sm.likelihood=localLikelihood;//+log(odo_gain);
			sm.pose=localPose;
			moveList.push_back(sm);
			//update the move list
		} while(move!=Done);
		currentPose=bestLocalPose;
		//cout << __PRETTY_FUNCTION__ << "currentScore=" << currentScore<< endl;
		//here we look for the best move;
	}while (currentScore>bestScore || refinement<m_optRecursiveIterations);
	//cout << __PRETTY_FUNCTION__ << "bestScore=" << bestScore<< endl;
	//cout << __PRETTY_FUNCTION__ << "iterations=" << count<< endl;
	
	//normalize the likelihood
	double lmin=1e9;
	double lmax=-1e9;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		lmin=it->likelihood<lmin?it->likelihood:lmin;
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	//cout << "lmin=" << lmin << " lmax=" << lmax<< endl;
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++){
		it->likelihood=exp(it->likelihood-lmax);
		//cout << "l=" << it->likelihood << endl;
	}
	//compute the mean
	OrientedPoint mean(0,0,0);
	double lacc=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		lacc+=it->likelihood;
	}
	mean=mean*(1./lacc);
	//OrientedPoint delta=mean-currentPose;
	//cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lacc, cov.xy/=lacc, cov.xt/=lacc, cov.yy/=lacc, cov.yt/=lacc, cov.tt/=lacc;
	
	_mean=currentPose;
	_cov=cov;
	return bestScore;
}

void ScanMatcher::setLaserParameters
	(unsigned int beams, double* angles, const OrientedPoint& lpose){
	/*if (m_laserAngles)
		delete [] m_laserAngles;
	*/
	assert(beams<LASER_MAXBEAMS);
	m_laserPose=lpose;
	m_laserBeams=beams;
	//m_laserAngles=new double[beams];
	memcpy(m_laserAngles, angles, sizeof(double)*m_laserBeams);	
}
	

double ScanMatcher::likelihood
	(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	ScoredMoveList moveList;
	
	for (double xx=-m_llsamplerange; xx<=m_llsamplerange; xx+=m_llsamplestep)
	for (double yy=-m_llsamplerange; yy<=m_llsamplerange; yy+=m_llsamplestep)
	for (double tt=-m_lasamplerange; tt<=m_lasamplerange; tt+=m_lasamplestep){
		
		OrientedPoint rp=p;
		rp.x+=xx;
		rp.y+=yy;
		rp.theta+=tt;
		
		ScoredMove sm;
		sm.pose=rp;
		
		likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
		moveList.push_back(sm);
	}
	
	//OrientedPoint delta=mean-currentPose;
	//cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
	//normalize the likelihood
	double lmax=-1e9;
	double lcum=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++){
		//it->likelihood=exp(it->likelihood-lmax);
		lcum+=exp(it->likelihood-lmax);
		it->likelihood=exp(it->likelihood-lmax);
		//cout << "l=" << it->likelihood << endl;
	}
	
	OrientedPoint mean(0,0,0);
	double s=0,c=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		s+=it->likelihood*sin(it->pose.theta);
		c+=it->likelihood*cos(it->pose.theta);
	}
	mean=mean*(1./lcum);
	s/=lcum;
	c/=lcum;
	mean.theta=atan2(s,c);
	
	
	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lcum, cov.xy/=lcum, cov.xt/=lcum, cov.yy/=lcum, cov.yt/=lcum, cov.tt/=lcum;
	
	_mean=mean;
	_cov=cov;
	_lmax=lmax;
	return log(lcum)+lmax;
}

double ScanMatcher::likelihood
	(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p,
	Gaussian3& odometry, const double* readings, double gain){
	ScoredMoveList moveList;
	
	
	for (double xx=-m_llsamplerange; xx<=m_llsamplerange; xx+=m_llsamplestep)
	for (double yy=-m_llsamplerange; yy<=m_llsamplerange; yy+=m_llsamplestep)
	for (double tt=-m_lasamplerange; tt<=m_lasamplerange; tt+=m_lasamplestep){
		
		OrientedPoint rp=p;
		rp.x+=xx;
		rp.y+=yy;
		rp.theta+=tt;
		
		ScoredMove sm;
		sm.pose=rp;
		
		likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
		sm.likelihood+=odometry.eval(rp)/gain;
		assert(!isnan(sm.likelihood));
		moveList.push_back(sm);
	}
	
	//OrientedPoint delta=mean-currentPose;
	//cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
	//normalize the likelihood
  double lmax=-std::numeric_limits<double>::max();
	double lcum=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++){
		//it->likelihood=exp(it->likelihood-lmax);
		lcum+=exp(it->likelihood-lmax);
		it->likelihood=exp(it->likelihood-lmax);
		//cout << "l=" << it->likelihood << endl;
	}
	
	OrientedPoint mean(0,0,0);
	double s=0,c=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		s+=it->likelihood*sin(it->pose.theta);
		c+=it->likelihood*cos(it->pose.theta);
	}
	mean=mean*(1./lcum);
	s/=lcum;
	c/=lcum;
	mean.theta=atan2(s,c);
	
	
	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lcum, cov.xy/=lcum, cov.xt/=lcum, cov.yy/=lcum, cov.yt/=lcum, cov.tt/=lcum;
	
	_mean=mean;
	_cov=cov;
	_lmax=lmax;
	double v=log(lcum)+lmax;
	assert(!isnan(v));
	return v;
}

void ScanMatcher::setMatchingParameters
	(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations,  double likelihoodSigma, unsigned int likelihoodSkip){	
	m_usableRange=urange;
	m_laserMaxRange=range;
	m_kernelSize=kernsize;
	m_optLinearDelta=lopt;
	m_optAngularDelta=aopt;
	m_optRecursiveIterations=iterations;
	m_gaussianSigma=sigma;
	m_likelihoodSigma=likelihoodSigma;
	m_likelihoodSkip=likelihoodSkip;
}

};

