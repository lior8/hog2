/*
 *  $Id: sample.cpp
 *  hog2
 *
 *  Created by Nathan Sturtevant on 5/31/05.
 *  Modified by Nathan Sturtevant on 02/29/20.
 *
 * This file is part of HOG2. See https://github.com/nathansttt/hog2 for licensing information.
 *
 */

#include "Common.h"
#include "Driver.h"
#include "UnitSimulation.h"
#include "EpisodicSimulation.h"
#include "Map2DEnvironment.h"
#include "RandomUnits.h"
#include "TemplateAStar.h"
#include "GraphEnvironment.h"
#include "ScenarioLoader.h"
#include "BFS.h"
#include "PEAStar.h"
#include "EPEAStar.h"
#include "MapGenerators.h"
#include "FPUtil.h"
#include "CanonicalGrid.h"

bool mouseTracking = false;
bool runningSearch1 = false;
bool startNewProblem = true;
int gStepsPerFrame = 1;
double searchWeight = 0;
bool setTarget = true;
bool recording = false;
void LoadMap(Map *m, int which);
void StartSearch();

TemplateAStar<xyLoc, tDirection, MapEnvironment> a1;
MapEnvironment *me = 0;

std::vector<xyLoc> path;
std::vector<xyLoc> targets;
xyLoc start;

int main(int argc, char* argv[])
{
	InstallHandlers();
	RunHOGGUI(argc, argv, 1000, 1000);
	return 0;
}


/**
 * Allows you to install any keyboard handlers needed for program interaction.
 */
void InstallHandlers()
{
//	InstallKeyboardHandler(MyDisplayHandler, "Cycle Abs. Display", "Cycle which group abstraction is drawn", kAnyModifier, '\t');
	InstallKeyboardHandler(MyDisplayHandler, "Record", "Record the screen.", kNoModifier, 'r');
	InstallKeyboardHandler(MyDisplayHandler, "Pause Simulation", "Pause simulation execution.", kNoModifier, 'p');
	InstallKeyboardHandler(MyDisplayHandler, "Step Simulation", "If the simulation is paused, step forward .1 sec.", kNoModifier, 'o');
	InstallKeyboardHandler(MyDisplayHandler, "Set target", "Add target locations", kNoModifier, 't');
	InstallKeyboardHandler(MyDisplayHandler, "Set start", "Add start location", kNoModifier, 's');
	InstallKeyboardHandler(MyDisplayHandler, "Dijkstra", "Use Dijkstra", kNoModifier, '0');
	InstallKeyboardHandler(MyDisplayHandler, "A*", "Optimal A*", kNoModifier, '1');
	InstallKeyboardHandler(MyDisplayHandler, "WA*(2)", "A* with weight of 2", kNoModifier, '2');
	InstallKeyboardHandler(MyDisplayHandler, "WA*(10)", "A* with weight of 10", kNoModifier, '3');
	InstallKeyboardHandler(MyDisplayHandler, "WA*(100)", "A* with weight of 100", kNoModifier, '4');
	InstallKeyboardHandler(MyDisplayHandler, "Forward Dijkstra", "Dijkstra in the forward direction", kNoModifier, '5');
	InstallKeyboardHandler(MyDisplayHandler, "Faster", "Increase simulation speed", kAnyModifier, ']');
	InstallKeyboardHandler(MyDisplayHandler, "Slower", "Decrease simulation speed", kAnyModifier, '[');
	
	InstallWindowHandler(MyWindowHandler);
	
	InstallMouseClickHandler(MyClickHandler);
}

void MyWindowHandler(unsigned long windowID, tWindowEventType eType)
{
	if (eType == kWindowDestroyed)
	{
		printf("Window %ld destroyed\n", windowID);
		RemoveFrameHandler(MyFrameHandler, windowID, 0);
		
		runningSearch1 = false;
		mouseTracking = false;
	}
	else if (eType == kWindowCreated)
	{
		printf("Window %ld created\n", windowID);
		InstallFrameHandler(MyFrameHandler, windowID, 0);
		ReinitViewports(windowID, {-1, -1, 1, 1}, kScaleToSquare);

		Map *map = new Map(1, 1);
		LoadMap(map, 0);
		map->SetTileSet(kWinter);
		me = new MapEnvironment(map);
		me->SetDiagonalCost(1.5);
	}
	
}

const char *GetAlgName()
{
	switch ((int)searchWeight)
	{
		case -1: return "Forward Dijkstra";
		case 0: return "Backwards Dijkstra";
		case 1: return "Backwards A*";
		case 2: return "Backwards WA*(2)";
		case 10: return "Backwards WA*(10)";
		case 100: return "Backwards WA*(100)";
		default: return "";
	}
}

void MyFrameHandler(unsigned long windowID, unsigned int viewport, void *)
{
	Graphics::Display &display = getCurrentContext()->display;
	std::string s;
	static bool init = false;
	if (init == false)
	{
		display.StartBackground();
		me->Draw(display);
		display.EndBackground();
		init = true;
	}

	if (mouseTracking)
	{
		// me->SetColor(Colors::red); 
		// me->DrawLine(display,xyLoc(px1, py1), xyLoc(px2, py2));
	}
	
	if (runningSearch1)
	{
		for (int x = 0; x < gStepsPerFrame; x++)
		{
			if (a1.DoSingleSearchStep(path))
			{
//					printf("Solution: moves %d, length %f, %lld nodes, %u on OPEN\n",
//						   (int)path.size(), ma1->GetPathLength(path), a1.GetNodesExpanded(),
//						   a1.GetNumOpenItems());
				s = GetAlgName();
				s += ": "+std::to_string(a1.GetNodesExpanded()) + " expansions";
				//snprintf(buffer, 255, "%" PRId64 " expansions", a1.GetNodesExpanded());
				submitTextToBuffer(s.c_str());
				runningSearch1 = false;
				break;
			}
			if (searchWeight == -1) // forward dijsktra
			{
				bool found = false;
				double cost;
				for (int y = 0; y < targets.size(); y++)
				{
					if (a1.GetClosedListGCost(targets[y], cost))
					{
						found = true;
						a1.ExtractPathToStart(targets[y], path);
					}
				}
				if (found)
				{
					s = GetAlgName();
					s += ": "+std::to_string(a1.GetNodesExpanded()) + " expansions";
//					snprintf(buffer, 255, "%" PRId64 " expansions", a1.GetNodesExpanded());
					submitTextToBuffer(s.c_str());
					runningSearch1 = false;
					break;
				}
			}
		}
	}
	else { // not running search
		for (auto i : targets)
			me->Draw(display, i);
	}
	// hack to handle that forward search doesn't
	
	if (startNewProblem == true)
		a1.Draw(display);

	if (path.size() != 0)
	{
		me->SetColor(0, 0.5, 1);
		for (int x = 1; x < path.size(); x++)
		{
			me->DrawLine(display, path[x-1], path[x], 10);
		}

	}
}



void MyDisplayHandler(unsigned long windowID, tKeyboardModifier mod, char key)
{
	char messageString[255];
	switch (key)
	{
		case 't':
			setTarget = true;
			break;
		case 's':
			if (targets.size() > 1)
			{
				submitTextToBuffer("Choose agent location.");
				setTarget = false;
			}
			else {
				submitTextToBuffer("Need to set two or more target locations first.");
			}
			break;
		case 'r':
			a1.SetReopenNodes(!a1.GetReopenNodes());
			if (a1.GetReopenNodes())
				submitTextToBuffer("Reopenings enabled");
			else
				submitTextToBuffer("Reopenings disabled");
			break;
		case '0':
			searchWeight = 0;
			submitTextToBuffer("Algorithm: Dijkstra");
			if (runningSearch1 || path.size() > 0)
				StartSearch();
			break;
		case '1':
			searchWeight = 1;
			submitTextToBuffer("Algorithm: A*");
			if (runningSearch1 || path.size() > 0)
				StartSearch();
			break;
		case '2':
			searchWeight = 2;
			submitTextToBuffer("Algorithm: WA*(2)");
			if (runningSearch1 || path.size() > 0)
				StartSearch();
			break;
		case '3':
			searchWeight = 10;
			submitTextToBuffer("Algorithm: WA*(10)");
			if (runningSearch1 || path.size() > 0)
				StartSearch();
			break;
		case '4':
			searchWeight = 100;
			submitTextToBuffer("Algorithm: WA*(100)");
			if (runningSearch1 || path.size() > 0)
				StartSearch();
			break;
		case '5':
			searchWeight = -1;
			submitTextToBuffer("Algorithm: Dijkstra Forward");
			if (runningSearch1 || path.size() > 0)
				StartSearch();
			break;
		case '[':
			gStepsPerFrame /= 2;
			snprintf(messageString, 255, "Speed: %d steps per frame", gStepsPerFrame);
			submitTextToBuffer(messageString);
			break;
		case ']':
			if (gStepsPerFrame <= 16384)
				gStepsPerFrame *= 2;
			if (gStepsPerFrame == 0)
				gStepsPerFrame = 1;
			snprintf(messageString, 255, "Speed: %d steps per frame", gStepsPerFrame);
			submitTextToBuffer(messageString);
			break;
//		case 'p': unitSims[windowID]->SetPaused(!unitSims[windowID]->GetPaused()); break;
		case 'o':
		{
		}
			break;
		default:
			break;
	}
}

void StartSearch()
{
	if (targets.size() == 0)
		return;
	a1.SetStopAfterGoal(true);
	a1.SetReopenNodes(false);
	me->SetEightConnected();

	if (searchWeight == -1)
	{
		a1.SetWeight(0);
		a1.InitializeSearch(me, start, targets[0], path);
	}
	else {
		a1.SetWeight(searchWeight);
		a1.InitializeSearch(me, targets[0], start, path);
		for (auto x = 1; x < targets.size(); x++)
			a1.AddAdditionalStartState(targets[x]);
	}
	mouseTracking = false;
	runningSearch1 = true;
}

bool MyClickHandler(unsigned long windowID, int, int, point3d loc, tButtonType button, tMouseEventType mType)
{
	if (setTarget == false)
	{
		if (mType == kMouseDown)
		{
			int x, y;
			me->GetMap()->GetPointFromCoordinate(loc, x, y);
			if ((me->GetMap()->GetTerrainType(x, y)>>terrainBits) == (kGround>>terrainBits))
			{
				start.x = (uint16_t)x;
				start.y = (uint16_t)y;
				StartSearch();
				startNewProblem = true;
			}
			return true;
		}
	}
	if (setTarget == true)
	{
		switch (mType)
		{
			case kMouseDown:
			{
				if (startNewProblem)
					targets.resize(0);
				startNewProblem = false;
				runningSearch1 = false;
				path.resize(0);
				int x, y;
				me->GetMap()->GetPointFromCoordinate(loc, x, y);
				if ((me->GetMap()->GetTerrainType(x, y)>>terrainBits) != (kGround>>terrainBits))
					break;
				xyLoc next(x, y);
				targets.push_back(next);
			}
			default:
				break;
		}
		return true;
	}
	return false;
}




void BGMap2(Map *m)
{
	m->Scale(148, 139);
	const char map[] = "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.....@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.......@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@..........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...@@@@..........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...@@@@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...@@@@....@@@.......@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...@@@@@@....@@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@@@@@@................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.....@@@@@@@@.@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.......@.@@@@..@@...............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@...@@@..@@@................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@....@..@@@.................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@....................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@.....................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@..@.......................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@@...........................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@........@@.................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....................@@@@...............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...................@@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.................@@@@@@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@@@.@@@@@@@...............@@@@@@@@..........@@@@@..@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@@.......@@@@@.............@@@@@@@@..........@@@@@...@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@.........@@@@@...........@@@@@@@@........@@@@@@@........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@@..........@@@@@@@@.........@@@@@@@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@..........@@@@@@@.........@@@@@@...........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@.........@@@@@@@.........@@@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@..........@@@@@@.........@@@@@@.....@@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@.............@..@.........@@@@@@.....@@@@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@..........................@@@@@@.....@@@@@..........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@.........................@@@@@@......@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@..........................@@@@@@.....@@@@...............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@..........................@@@@@@.....@@@@.................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@..........................@@@@@@.....@@@@...................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........................@@.................@@@@@@.....@@@@...@@@...............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.............@@@.........@@@@...............@@@@@@@....@@@@...@@@@@......@........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............@@@@@.......@@@@@@.............@@@@@......@@@@...@@@@@@.@...@@@........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............@@@@@.......@@@@@@@............@@@@@......@@@@...@@@@@@@@@..@@@@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..........@@@@@.......@@@@@@@@............@@@@......@@@@....@@@@@@@@..@@@@...........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.........@@@@@.......@@@@@@@@@.............@@@.......@@.....@@@@@@@..@@@@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.......@@@@@.......@@@@@@@@................@@................@@.@..@@@@...............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.....@@@@@.......@@@@@@@@..........@.............................@@@@..................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...@@@@@........@@@@@@@..........@@@...........................@@@@..@@...............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...@@@@@@@..@@@@.........@@@@@@..........@@@@@.........................@@@@..@@@@......@........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@@@@@@@@@..........@@@@@..........@@@@@@........................@@@@..@@@@@@....@@@........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@......@@@@@@@@............@@@..........@@@@@@@.......................@@@@..@@@@@@@...@@@@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@......@@@@@@.........................@@@@@@........................@@@@..@@@@@@@...@@@@...........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.......@@@@@@@......................@@@@@@....................@@...@@@....@@@@@...@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.........@@@@@@.....................@@@@@@....................@@@@...@......@@@...@@@@................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@.....@.....@@@@@@...................@@@@@@....................@@@@@@..............@@@@..................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@.....@@@.....@@@@@..................@@@@@@@........@.@@.......@@@..@@@............@@@@..@@@................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.........@@@@@.....@@@@.................@@@@@..........@@@@@@.....@@@....@@@..........@@@@@@@@@@@................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........@@@@@@.....@..................@@@@@...........@@@@@@....@@@......@@@........@@@@@@@@@@@@.......@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........@@@@@@........@@@...........@@@@@...........@@@@@@....@@@........@@@......@@@@.@@@@@@@@@.....@@@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........@@@@@........@@@@.........@@@@@@............@@@@....@@@..........@@@....@@@@..@.@@@@@@@....@@@@..........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........@@@..@......@@@@@.......@@@@@@@.............@@....@@@............@@@...@@@.....@@@@@@....@@@@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........@..@@@.......@@@@.....@@@@@@@@@.................@@@..............@@....@.......@@@@....@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.........@@@@@.......@@@@...@@@@@@@@...................@@...............@@@............@@....@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........@@@@@@.......@@@@.@@@@@@@@...................@@.........@@......@@@................@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.........@@@@.......@@@@@@@@@@@@...................@@@.........@@@.....@@@...............@@@@.............@@@@@@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.........@@.......@@@@@@@@@@@@......@@...........@@@@..........@@......@...............@@@@.............@@@@@@....@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@................@@@@@@@@@@@@......@@@@.........@@@@@@.........@@.....................@@@@.............@@@@@@@...@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@@@@@@@@@@......@@@@@@........@@@@@@@........@@@.....@@............@@@@.............@@@@@@@@...@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............@@@@@@@@@@@@@.....@@@@@@@@@......@@@@@@@@@@....@@@@.....@@...........@@@@.............@@@@@@@@@....@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..........@@@@@@@@@@@@@@....@@@@@@@@@@@......@@@@@@@@@@@@@@..@@.................@@@.............@@@@@@@@.......@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........@@@@@@@@@@@@....@@@@@@@@@@@@@@@S.....@@@@@@@@@@.....@@..................@.............@@@@@@@@.........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@......@@@@@@@@@@@@.....@@@@@@@@@@@@@@@SS.....@@@@@@@@.......@...............................@@@@@@@@@..........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@@@@@@@@@@......@@@@@@@@@@@@@@SSS.......@@@@@.......................................@@@@@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@@@@@@.@@........@@@@@@@@@@@@SSSS........@@@.......................................@@@@@@@@...............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............@@@@@@@@@@SSSSS@@.......@.......................................@@@@@@@@.................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............@@@@@@@@@SSSSS@@@@@............................................@@@@@@@@...................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@@@SSSSSS@@@@@@@..........................................@@@@@@@@.....................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@@@@SSSSS@@@@@@@@@........................................@@@@@@@@..........@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@@@@@SSSS@@@@@@@@@@@......................................@@@@@@@@..........@@@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@@@@@@SSSS@@@@@@@@@@@@S...................................@@@@@@@@..........@@@@@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.............@@@@@@@..SSS@@@@@@@@@@@@@@@SS................................@@@@@@@@............@@@@...........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@@@@@...SS@@@@@@@@@@@@@@@SSS...............................@@@@@@@@..............@@...........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.........@@@@@@@@.....S@@@@@@@@@@@@SSSS................................@@@@@@@...........................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@......@@@@@@...........@@@@@@@@@@SSSSS@@...............................@@.@@@.....@....................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@@@@@...........@@@@@@@@@SSSSS@@@@...............................@........@@@..................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@@@@.............@@@@@@@SSSSS@@@@@@......................................@@@@@................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@@@@SSSSS@@@@@@@@......................................@@@@@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@@@@SSSS@@@@@@@@@@......................................@@@@@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@@@@@SSS@@@@@@@@@@@@......................................@@@@@..........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@@@@SSSS@@@@@@@@@@@@@@......................................@@@@@........@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............@@@@@@@......@@@@@@@@@@@@@@................................@......@@@@@......@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..........@@@@@@@@.....@@@@@@@@@@@@@@.............@...................@@......@@@@@....@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.......@@@@@@.@@@......@@@@@@@@@@@@.............@@@...................@@......@@@@@..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.....@@@@@@@...........@@@@@@@@@@.............@@@@@...................@@......@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...@@@@@@@@...........@@@@@@@@@.............@@@@@@...........@@.......@.......@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.@@@@@@@.............@@@@@@@@.............@@@@@@@..........@@@@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.............@@@@@@@@.............@@@@@@............@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.............@@@@@@@@@............@@@@@@..............@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............@@@@@@@@@@...........@@@@@@@............................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.............@@@@@@@@.....@.......@@@@@@@@...........................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@@@@@@......@@.....@@@@@@@@...........@@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.........@@@@@@@@.......@@@...@@@@@@@@.............@@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.......@@@@@@...........@@@@@@@@@@@@...............@@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@@@@@@...........@@@@@@@@@@@.......@@........@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@@@@@............@@@@@@@@@@.......@@@@........@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@@@@@@@@@........@@@@...................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@@@@@@@............@@@@.................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@@@@@@@..............@@@@...............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...........@@@@@@@@@................@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............@@@@@@@@...................@@@@......@@@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..........@@@@@@@@.....................@@@@.....@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........@@@@@@@@.......................@@@@...@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@......@@@@@@@@.........................@@@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@@@@@@@..........................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...@@@@@.@@...........@@@..............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.@@@@@..............@@@@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@@@@............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.................@@.............@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@................................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..............................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@............................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..........................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@........................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@......................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....................@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
	int which = 0;
	for (int y = 0; y < m->GetMapHeight(); y++)
	{
		for (int x = 0; x < m->GetMapWidth(); x++)
		{
			if (map[which] == '.')
				m->SetTerrainType(x, y, kGround);
			else if (map[which] == 'S')
				m->SetTerrainType(x, y, kSwamp);
			else
				m->SetTerrainType(x, y, kOutOfBounds);
			which++;
		}
	}
//	mapChange = true;
}


void LoadMap(Map *m, int which)
{
	assert(m != 0);
	BGMap2(m);
}

