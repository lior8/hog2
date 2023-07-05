#include <cstring>
#include "PermutationPDB.h"
#include "LexPermutationPDB.h"
#include "MR1PermutationPDB.h"
#include "MNPuzzle.h"
#include "TOH.h"
#include "Timer.h"
#include "STPInstances.h"
#include "TemplateAStar.h"
#include "TASA.h"
#include "TAS.h"
#include <iostream>
#include <vector>
#include "MultipleAdditiveHeuristic.h"
#include "TASS.h"
#include "FastDict.h"

using namespace std;


template<class state>
class HeuristicsGroup
{
public:
	mutable vector<Heuristic<state>*> heuristics;
	mutable int count;
	HeuristicsGroup(vector<Heuristic<state>*> heu, int _count)
	{
		count = _count;
		heuristics = heu;
	}
	HeuristicsGroup()
	{
		heuristics.clear();
		count = 0;
	}
	~HeuristicsGroup()
	{

	}
};


template<class state>
class GroupHeuristic : public Heuristic<state>
{
    public:
	vector<HeuristicsGroup<state>*> groups;
    double HCost(const state &state1, const state &state2) const
    {
        double res = 0;
		int bestIndex = 0;
        for (int j = 0; j < groups.size(); j++)
        {
			auto group = groups[j];
            double partialRes = 0;
            for (auto h : group->heuristics)
            {
                partialRes += fabs(h->HCost(state1, state1) - h->HCost(state2, state2));
            }
			//res += partialRes;
            res = max(res, partialRes);
			if (res == partialRes)
			{
				bestIndex = j;
				//groups[j].second++;
			}
        }
		groups[bestIndex]->count++;// = HeuristicGroup(groups[bestIndex].indices, groups[bestIndex].count);
	    return res;
    	//return abs(baseH.HCost(state1, state1) - baseH.HCost(state2, state2));
    }

    void AddGroup(vector<Heuristic<state>*> _heuristics)
    {
		HeuristicsGroup<state> *hg = new HeuristicsGroup<state>(_heuristics, 0);
        groups.push_back(hg);
    }

	void RemoveGroup(int index)
	{
		groups.erase(groups.begin() + index);
	}

	void ClearCounts()
	{
		for (int i = 0; i < groups.size(); i++)
		{
			groups[i]->count = 0;
		}
	}

	int GetLeastAccessed()
	{
		int min = 1000000000;
		int ind = -1;
		for (int i = 0; i < groups.size(); i++)
		{
			//cout << "count " << i << ": " << groups[i]->count << endl;
			if (groups[i]->count < min)
			{
				min = groups[i]->count;
				ind = i;
			}
		}
		//cout << "min: " << ind << endl;
		return ind;
	}
};


template <int numDisks, int pdb1Disks, int pdb2Disks = numDisks-pdb1Disks>
Heuristic<TOHState<numDisks>> *BuildPDB(const TOHState<numDisks> &goal)
{
	TOH<numDisks> toh;
	TOH<pdb1Disks> absToh1;
	TOH<pdb2Disks> absToh2;
	TOHState<pdb1Disks> absTohState1;
	TOHState<pdb2Disks> absTohState2;
	
	
	TOHPDB<pdb1Disks, numDisks, pdb2Disks> *pdb1 = new TOHPDB<pdb1Disks, numDisks, pdb2Disks>(&absToh1, goal); // top disks
	TOHPDB<pdb2Disks, numDisks> *pdb2 = new TOHPDB<pdb2Disks, numDisks>(&absToh2, goal); // bottom disks
	pdb1->BuildPDB(goal, std::thread::hardware_concurrency(), false);
	pdb2->BuildPDB(goal, std::thread::hardware_concurrency(), false);
	
	Heuristic<TOHState<numDisks>> *h = new Heuristic<TOHState<numDisks>>;
	
	h->lookups.resize(0);
	
	h->lookups.push_back({kAddNode, 1, 2});
	h->lookups.push_back({kLeafNode, 0, 0});
	h->lookups.push_back({kLeafNode, 1, 1});
	h->heuristics.resize(0);
	h->heuristics.push_back(pdb1);
	h->heuristics.push_back(pdb2);
	
	return h;
}

template <int numDisks, int pdb1Disks, int pdb2Disks = numDisks-pdb1Disks>
void BuildSinglePDB(Heuristic<TOHState<numDisks>> &h, const TOHState<numDisks> &pivot1)
{
	TOH<numDisks> toh;
	TOH<pdb1Disks> absToh1;
	TOH<pdb2Disks> absToh2;
	TOHState<pdb1Disks> absTohState1;
	TOHState<pdb2Disks> absTohState2;
	
	
	//TOHPDB<pdb1Disks, numDisks, pdb2Disks> *pdb1 = new TOHPDB<pdb1Disks, numDisks, pdb2Disks>(&absToh1, pivot1); // top disks
	TOHPDB<pdb1Disks, numDisks> *pdb2 = new TOHPDB<pdb1Disks, numDisks>(&absToh1, pivot1); // bottom disks
	pdb2->BuildPDB(pivot1, std::thread::hardware_concurrency(), false);
	h.lookups.resize(0);
	h.heuristics.resize(0);
	h.heuristics.push_back(pdb2);
}

template<int numDisks, int pdb1Disks, int maxPivots = 5>
void HUpdate(TOHState<numDisks> s, TOHState<numDisks> g, Heuristic<TOHState<numDisks>>*& h)
{
	GroupHeuristic<TOHState<numDisks>>* mhh;
	mhh = (GroupHeuristic<TOHState<numDisks>>*)h;
	//cout << mhh->groups.size() << endl;
	TOH<numDisks> toh;
	TOH<pdb1Disks> absToh1;
	TOHState<pdb1Disks> absTohState1;

	TOHPDB<pdb1Disks, numDisks> *pdb2 = new TOHPDB<pdb1Disks, numDisks>(&absToh1, s); // bottom disks
	TOHPDB<pdb1Disks, numDisks> *pdb4 = new TOHPDB<pdb1Disks, numDisks>(&absToh1, g); // bottom disks
	pdb2->BuildPDB(s, std::thread::hardware_concurrency(), false);
	pdb4->BuildPDB(g, std::thread::hardware_concurrency(), false);

	Heuristic<TOHState<numDisks>> __h;
	while (mhh->groups.size() + 2 > maxPivots)
	{
		auto la = mhh->GetLeastAccessed();
		mhh->RemoveGroup(la);
	}
	mhh->AddGroup({pdb2});
	mhh->AddGroup({pdb4});
	mhh->ClearCounts();
	h = mhh;
}


double Phi(double h, double g)
{
    return h;
}

void Test5(int problemSeed)
{
    srandom(problemSeed);
	//TemplateAStar<MNPuzzleState<5, 5>, slideDir, MNPuzzle<5,5>> astar;
	MNPuzzle<5,5> mnp;
	MNPuzzleState<5, 5> start, goal;
    vector<MNPuzzleState<5, 5>> path1, path2;
    start = mnp.Generate_Random_Puzzle();
    goal = mnp.Generate_Random_Puzzle();
    cout << mnp.HCost(start, goal) << endl;
    cout << start << endl;
    cout << goal << endl;

    TemplateAStar<MNPuzzleState<5, 5>, slideDir, MNPuzzle<5,5>> astar;
    astar.SetPhi(Phi);
    TASA<MNPuzzle<5,5>, MNPuzzleState<5,5>> tasa(&mnp, start, goal, &mnp, &mnp, 10);
    Timer timer1, timer2;
    
    timer1.StartTimer();
    astar.GetPath(&mnp, start, goal, path1);
    timer1.EndTimer();
    cout << "GBFS: " << astar.GetNodesExpanded() << " " << timer1.GetElapsedTime() << " " << path1.size() << endl;

    timer2.StartTimer();
    tasa.GetPath(path2);
    timer2.EndTimer();
    cout << "TASA: " << tasa.GetNodesExpanded() << " " << timer2.GetElapsedTime() << " " << path2.size() << endl;
}


template<int disks>
void TOHTest(int problemSeed)
{
    srandom(problemSeed);
	//TemplateAStar<MNPuzzleState<5, 5>, slideDir, MNPuzzle<5,5>> astar;
	TOH<disks> toh;
	TOHState<disks> start, goal;
    vector<TOHState<disks>> path1, path2;
    MultipleAdditiveHeuristic<TOHState<disks>> h;

    start.counts[0] = start.counts[1] = start.counts[2] = start.counts[3] = 0;
	for (int x = disks; x > 0; x--)
	{
		int whichPeg = random()%4;
		start.disks[whichPeg][start.counts[whichPeg]] = x;
		start.counts[whichPeg]++;
	}
	goal.counts[0] = goal.counts[1] = goal.counts[2] = goal.counts[3] = 0;
	for (int x = disks; x > 0; x--)
	{
		int whichPeg = random()%4;
		goal.disks[whichPeg][goal.counts[whichPeg]] = x;
		goal.counts[whichPeg]++;
	}

    cout << start << endl;
    cout << goal << endl;

    auto fh = BuildPDB<disks, disks / 2>(goal);
    auto bh = BuildPDB<disks, disks / 2>(start);
    h.AddHeuristic(fh->heuristics[0]);
    h.AddHeuristic(fh->heuristics[1]);
    h.AddHeuristic(bh->heuristics[0]);
    h.AddHeuristic(bh->heuristics[1]);
    h.AddGroup({0, 1});
    h.AddGroup({2, 3});


    TemplateAStar<TOHState<disks>, TOHMove, TOH<disks>> astar;
    astar.SetPhi(Phi);
    astar.SetHeuristic(&h);
    TASA<TOH<disks>, TOHState<disks>> tasa(&toh, start, goal, &h, &h, 10);
    Timer timer1, timer2;
    
    timer1.StartTimer();
    astar.GetPath(&toh, start, goal, path1);
    timer1.EndTimer();
    cout << "GBFS: " << astar.GetNodesExpanded() << " " << timer1.GetElapsedTime() << " " << path1.size() << endl;

    timer2.StartTimer();
    tasa.GetPath(path2);
    timer2.EndTimer();
    cout << "TASA: " << tasa.GetNodesExpanded() << " " << timer2.GetElapsedTime() << " " << path2.size() << endl;
}

/*
template<int disks>
void TASOptTest(int problemSeed)
{
    srandom(problemSeed);
	//TemplateAStar<MNPuzzleState<5, 5>, slideDir, MNPuzzle<5,5>> astar;
	TOH<disks> toh;
	TOHState<disks> start, goal;
    vector<TOHState<disks>> path1, path2;
    MultipleAdditiveHeuristic<TOHState<disks>> h;

    start.counts[0] = start.counts[1] = start.counts[2] = start.counts[3] = 0;
	for (int x = disks; x > 0; x--)
	{
		int whichPeg = random()%4;
		start.disks[whichPeg][start.counts[whichPeg]] = x;
		start.counts[whichPeg]++;
	}
	goal.counts[0] = goal.counts[1] = goal.counts[2] = goal.counts[3] = 0;
	for (int x = disks; x > 0; x--)
	{
		int whichPeg = random()%4;
		goal.disks[whichPeg][goal.counts[whichPeg]] = x;
		goal.counts[whichPeg]++;
	}

    cout << start << endl;
    cout << goal << endl;

    auto fh = BuildPDB<disks, disks / 2>(goal);
    auto bh = BuildPDB<disks, disks / 2>(start);
    h.AddHeuristic(fh->heuristics[0]);
    h.AddHeuristic(fh->heuristics[1]);
    h.AddHeuristic(bh->heuristics[0]);
    h.AddHeuristic(bh->heuristics[1]);
    h.AddGroup({0, 1});
    h.AddGroup({2, 3});


    TemplateAStar<TOHState<disks>, TOHMove, TOH<disks>> astar;
    astar.SetPhi(Phi);
    astar.SetHeuristic(&h);
    //TASA<TOH<disks>, TOHState<disks>> tasa(&toh, start, goal, &h, &h, 10);
    TASOpt<TOH<disks>, TOHState<disks>> tasOpt(&toh, start, goal, &h, &h, 10);
    Timer timer1, timer2;

    //timer2.StartTimer();
    //tasa.GetPath(path2);
    //timer2.EndTimer();
    //cout << "TASA: " << tasa.GetNodesExpanded() << " " << timer2.GetElapsedTime() << " " << path2.size() << endl;
    

    timer1.StartTimer();
    astar.GetPath(&toh, start, goal, path1);
    timer1.EndTimer();
    cout << "GBFS: " << astar.GetNodesExpanded() << " " << timer1.GetElapsedTime() << " " << path1.size() << endl;

    timer2.StartTimer();
    tasOpt.GetPath(path2);
    timer2.EndTimer();
    cout << "TASO: " << tasOpt.GetNodesExpanded() << " " << timer2.GetElapsedTime() << " " << path2.size() << endl;
    

    delete fh;
    delete bh;
}
*/

template<int numOfDisks, int pdb1Disks>
void PlanningTest(int seed_problem, int samples, kAnchorSelection selection, int episode, GroupHeuristic<TOHState<numOfDisks>> hm)
{
	std::cout << "TAS(" << samples << ")" << std::endl;
    srandom(seed_problem);
    TOH<numOfDisks> *toh = new TOH<numOfDisks>();
    TOHState<numOfDisks> start, goal, pivot1, pivot2;
    Heuristic<TOHState<numOfDisks>> hf, hb, h1, h2, h3, h4;
    vector<TOHState<numOfDisks>> path, aPath, fPath, bPath;

    start.counts[0] = start.counts[1] = start.counts[2] = start.counts[3] = 0;
	for (int x = numOfDisks; x > 0; x--)
	{
		int whichPeg = random()%4;
		start.disks[whichPeg][start.counts[whichPeg]] = x;
		start.counts[whichPeg]++;
	}
    goal.counts[0] = goal.counts[1] = goal.counts[2] = goal.counts[3] = 0;
	for (int x = numOfDisks; x > 0; x--)
	{
		int whichPeg = random()%4;
		goal.disks[whichPeg][goal.counts[whichPeg]] = x;
		goal.counts[whichPeg]++;
	}
	Heuristic<TOHState<numOfDisks>> spec1, spec2;
	BuildSinglePDB<numOfDisks, pdb1Disks>(spec1, start);
	BuildSinglePDB<numOfDisks, pdb1Disks>(spec2, goal);
	hm.AddGroup({spec1.heuristics[0]});
	hm.AddGroup({spec2.heuristics[0]});

    cout << start << endl;
    cout << goal << endl;
    std::cout << "------------------------------------" << std::endl;
	TASS<TOH<numOfDisks>, TOHState<numOfDisks>> tas(toh, start, goal, &hm, &hm, samples);
	tas.episode = episode;
	tas.SetHeuristicUpdate(HUpdate<numOfDisks, pdb1Disks, 2>);
    tas.SetAnchorSelection(selection);
	//astar.SetPhi(phi);
	//astar.SetHeuristic(&hm);
    Timer timer, timer2, timer3;
    timer.StartTimer();
	tas.GetPath(path);
    timer.EndTimer();
	cout << "Info:" << endl;
	cout << "Anchor Selection: " << selection << endl;
	cout << "Heuristic: " << endl;
	cout << "PDB size: " << pdb1Disks << endl;
	cout << "episode length: " << episode << endl;
	cout << "TAS: " << tas.GetNodesExpanded() << " " << timer.GetElapsedTime() << " " << path.size() << endl;
}




template<int numOfDisks, int pdb1Disks>
void GBFSTest(int seed_problem, GroupHeuristic<TOHState<numOfDisks>> hm)
{
    srandom(seed_problem);
    TOH<numOfDisks> *toh = new TOH<numOfDisks>();
    TOHState<numOfDisks> start, goal, pivot1, pivot2;
    Heuristic<TOHState<numOfDisks>> hf, hb, h1, h2, h3, h4;
    vector<TOHState<numOfDisks>> path, aPath, fPath, bPath;

    start.counts[0] = start.counts[1] = start.counts[2] = start.counts[3] = 0;
	for (int x = numOfDisks; x > 0; x--)
	{
		int whichPeg = random()%4;
		start.disks[whichPeg][start.counts[whichPeg]] = x;
		start.counts[whichPeg]++;
	}
    goal.counts[0] = goal.counts[1] = goal.counts[2] = goal.counts[3] = 0;
	for (int x = numOfDisks; x > 0; x--)
	{
		int whichPeg = random()%4;
		goal.disks[whichPeg][goal.counts[whichPeg]] = x;
		goal.counts[whichPeg]++;
	}
	Heuristic<TOHState<numOfDisks>> spec1, spec2;
	BuildSinglePDB<numOfDisks, pdb1Disks>(spec1, start);
	BuildSinglePDB<numOfDisks, pdb1Disks>(spec2, goal);
	hm.AddGroup({spec1.heuristics[0]});
	hm.AddGroup({spec2.heuristics[0]});

    cout << start << endl;
    cout << goal << endl;
    std::cout << "------------------------------------" << std::endl;
	TemplateAStar<TOHState<numOfDisks>, TOHMove, TOH<numOfDisks>> astar;
	astar.SetPhi(Phi);
	astar.SetHeuristic(&hm);
    Timer timer, timer2, timer3;
    timer.StartTimer();
	astar.GetPath(toh, start, goal, path);
    timer.EndTimer();
	cout << "Heuristic: " << endl;
	cout << "PDB size: " << pdb1Disks << endl;
	cout << "GBFS: " << astar.GetNodesExpanded() << " " << timer.GetElapsedTime() << " " << path.size() << endl;
}


void STPTest(int problemSeed)
{
    srandom(problemSeed);
	//TemplateAStar<MNPuzzleState<5, 5>, slideDir, MNPuzzle<5,5>> astar;
	MNPuzzle<4, 4> mnp;
	MNPuzzleState<4, 4> start, goal;
    vector<MNPuzzleState<4, 4>> path1, path2;
    start = mnp.Generate_Random_Puzzle();
    goal = mnp.Generate_Random_Puzzle();
    cout << mnp.HCost(start, goal) << endl;
    cout << start << endl;
    cout << goal << endl;

    TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
    astar.SetPhi(Phi);
    TASA<MNPuzzle<4, 4>, MNPuzzleState<4, 4>> tasa(&mnp, start, goal, &mnp, &mnp, 10);
    Timer timer1, timer2;
    
    timer1.StartTimer();
    astar.GetPath(&mnp, start, goal, path1);
    timer1.EndTimer();
    cout << "GBFS: " << astar.GetNodesExpanded() << " " << timer1.GetElapsedTime() << " " << path1.size() << endl;

    timer2.StartTimer();
    tasa.GetPath(path2);
    timer2.EndTimer();
    cout << "TASA: " << tasa.GetNodesExpanded() << " " << timer2.GetElapsedTime() << " " << path2.size() << endl;
}


int main(int argc, char* argv[])
{
	//TOHTest<18>(stoi(argv[1]));
    /*
    int size = stoi(argv[1]);
    if (size == 14)
    {
        for (int i = 1; i <= 100; i++)
            TASOptTest<14>(i);
    }
    else if (size == 16)
    {
        for (int i = 1; i <= 100; i++)
            TASOptTest<16>(i);
    }
    else if (size == 18)
    {
        TASOptTest<18>(stoi(argv[2]));
    }
    else if (size == 20)
    {
        TASOptTest<20>(stoi(argv[2]));
    }
    else if (size == 22)
    {
        TASOptTest<22>(stoi(argv[2]));
    }
    */
   	
   	int problemSize = stoi(argv[1]);
	if (problemSize == 20)
	{
		GroupHeuristic<TOHState<20>> gh;
    	PlanningTest<20, 12>(stoi(argv[2]), 10, (kAnchorSelection)stoi(argv[3]), stoi(argv[4]), gh);
	}
	else if (problemSize == 22)
	{
		GroupHeuristic<TOHState<22>> gh;
    	PlanningTest<22, 12>(stoi(argv[2]), 10, (kAnchorSelection)stoi(argv[3]), stoi(argv[4]), gh);
	}
	else if (problemSize == 24)
	{
		GroupHeuristic<TOHState<24>> gh;
    	PlanningTest<24, 12>(stoi(argv[2]), 10, (kAnchorSelection)stoi(argv[3]), stoi(argv[4]), gh);
	}
	else if (problemSize == 26)
	{
		GroupHeuristic<TOHState<26>> gh;
    	PlanningTest<26, 12>(stoi(argv[2]), 10, (kAnchorSelection)stoi(argv[3]), stoi(argv[4]), gh);
	}
	else if (problemSize == 28)
	{
		GroupHeuristic<TOHState<28>> gh;
    	PlanningTest<28, 12>(stoi(argv[2]), 10, (kAnchorSelection)stoi(argv[3]), stoi(argv[4]), gh);
	}
	else if (problemSize == 30)
	{
		GroupHeuristic<TOHState<30>> gh;
    	PlanningTest<30, 12>(stoi(argv[2]), 10, (kAnchorSelection)stoi(argv[3]), stoi(argv[4]), gh);
	}
	else if (problemSize == 32)
	{
		GroupHeuristic<TOHState<32>> gh;
    	PlanningTest<32, 12>(stoi(argv[2]), 10, (kAnchorSelection)stoi(argv[3]), stoi(argv[4]), gh);
	}
    
   	/*
   	GroupHeuristic<TOHState<24>> gh;
   	GBFSTest<24, 12>(stoi(argv[1]), gh);
	*/
	return 0;
}