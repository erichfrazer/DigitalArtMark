// Painting.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <windows.h>
#include <assert.h>
using namespace std;

struct StartStop
{
    int start;
    mutable int stop;

    StartStop()
    {
        start = 0;
        stop = 0;
    }

    StartStop(int s, int e)
    {
        start = s;
        stop = e;
    }

    inline bool operator < (const StartStop& other) const
    {
        return start < other.start;
    }
};

struct YandXStop
{
    int Y;
    mutable int XStop;

    YandXStop(int y, int xs)
    {
        Y = y;
        XStop = xs;
    }

    inline bool operator < (const YandXStop& other) const
    {
        return Y < other.Y;
    }
};

struct Y_HStroke
{
    int Y;
    mutable int XStart;
    mutable int XEnd;

    Y_HStroke()
    {

    }

    Y_HStroke(int y, int s, int e)
    {
        Y = y;
    }

    inline bool operator < (const Y_HStroke& other) const
    {
        return Y < other.Y;
    }

    inline bool operator == (Y_HStroke const& lhs)
    {
        // only one line in HLive can have the same Y
        return (lhs.Y == Y);
    }
};

struct X_YArray
{
    int X;
    vector<StartStop> Strokes;

    X_YArray()
    {
        X = 0;
    }

    ~X_YArray()
    {
    }

    X_YArray(int x, int n)
    {
        X = x;
        Strokes = vector<StartStop>(n);
    }

    inline bool operator < (const X_YArray& other) const
    {
        return X < other.X;
    }
};

typedef set<StartStop, less<StartStop>> StartStopSet;
typedef set<YandXStop, less<YandXStop>> YandXStopSet;
typedef map<int, StartStop> YHStrokeSet; // unordered is sig faster than ordered

map<int, StartStopSet> VStrokesRaw;     // used just to create VStrokes. Deleted after VStrokes is created
vector<X_YArray> VStrokes;              // array of X followed by Y start/stop pairs, sorted, no overlaps.
map<int, YandXStopSet> HStrokes;        // sorted by X, but possible overlaps
YHStrokeSet HLive;                      // lines that are CURRENT at SweepX, map key is Y value, StartStop are in X. Adjusted as we go.

list<pair<int, int>> foundPlusses;
int xmin = 0;
int xmax = 0;
int maxLinesInHLive = 0;

long binarySearchMap2(const vector<StartStop>& v, long L, long R, long long k)
{
    long M = (L + R) / 2;

    if (L > R)
    {
        return -1;
    }

    // found within the range
    if (k >= v[M].start && k <= v[M].stop)
    {
        return M;
    }

    if (v[M].start > k) // k is to the left of M
    {
        return binarySearchMap2(v, L, M - 1, k);
    }
    else // k is greater than M
    {
        return binarySearchMap2(v, M + 1, R, k);
    }
}

// this is ONLY called when we move on the xiterator for HStrokes. Should add linear time,
// since size of HLive is relatively constant, and # of lines at hstrokes_iter is relatively constant

void AddNewLinesToHLive(const map<int, YandXStopSet>::iterator& hstrokes_iter)
{
    int x = hstrokes_iter->first;

    //---------------------------------------------------------------------
    // add new lines to the Live set
    //---------------------------------------------------------------------

    const YandXStopSet& hs = hstrokes_iter->second;
    // there are a bunch of hstrokes in hs...all at different Y's
    // Can we find an existing line in Live and extend it?
    // (this means we need to search all existing live lines, no matter what their start x)
    for (const YandXStop& yxs : hs)
    {
        // This is messy having to do this. hs is sorted by y, so if we find the interator for HLive,
        // we should be able to increase it in lock-step and not have to call find() each time.
        // That would save a lot.

        YHStrokeSet::iterator i = HLive.find(yxs.Y);
        if (i == HLive.end()) // didnt find it, add it
        {
            HLive.emplace(pair<int, StartStop>(x, StartStop(x, yxs.XStop)));
        }
        else
        {
            // we need to insert a line at the same Y. If it's shorter, ignore it.
            // if it's longer, then extend the current one.

            if (x < i->second.start) continue;
            if (x > i->second.stop) continue;
            if (yxs.XStop <= i->second.stop) continue;
            i->second.stop = yxs.XStop;
        }
    }

    maxLinesInHLive = max(maxLinesInHLive, HLive.size());
    
}

void RemoveOldHLiveStrokes(int x)
{
    int minstopx = -1;

    YHStrokeSet::iterator hstroke_iter = HLive.begin();
    while (hstroke_iter != HLive.end())
    {
        if (hstroke_iter->second.stop < x)
        {
            hstroke_iter = HLive.erase(hstroke_iter);
        }
        else
        {
            if (minstopx == -1 || hstroke_iter->second.stop < minstopx)
            {
                minstopx = hstroke_iter->second.stop;
            }
            hstroke_iter++;
        }
    }
}

void CheckForPlusses(vector<X_YArray>::iterator vstrokeiter)
{
    int vertStrokeX = vstrokeiter->X;

    // this nugget's got to be FAST. Whatever you can do to speed this section up,
    // make it so.
    // VStrokeIter's YStroke list is sorted in Y
    // Likewise, the HLive iterator is ALSO sorted in Y!
    // Guess what? We can do this in nearly linear time!

    vector<StartStop>::iterator vstrokes_end = vstrokeiter->Strokes.end();
    vector<StartStop>::iterator vstrokes = vstrokeiter->Strokes.begin(); // ordered by Y
    YHStrokeSet::iterator hstroke_iter = HLive.begin(); // ordered by Y
    YHStrokeSet::iterator hstroke_iter_end = HLive.end();
    int hStrokesY;

    while (vstrokes != vstrokes_end)
    {
        // enumerate through hstrokes that are within range
        while (hStrokesY = hstroke_iter->first <= vstrokes->start)
        {
            hstroke_iter++;
            if (hstroke_iter == hstroke_iter_end)
            {
                return;
            }
        }

        // by there, hstrokesY should be "within range" of the vstroke. OR, it's way over.
        while (hStrokesY = hstroke_iter->first < vstrokes->stop)
        {
            // this should be forming a plus. We shouldn't have to do this test, we already know the horizontal line
            // crosses the vertical and makes a plus. But let's do it anyhow to make sure.
            if (hStrokesY > vstrokes->start && hStrokesY < vstrokes->stop)
            {
                foundPlusses.push_back(pair<int, int>(vertStrokeX, hStrokesY));
            }
            else
            {
                assert(false);
            }
            hstroke_iter++;
            if (hstroke_iter == hstroke_iter_end)
            {
                return;
            }
        }

        vstrokes++;
        if (vstrokes == vstrokes_end)
        {
            break;
        }
    }
}

long long getPlusSignCount(int N, vector<int> L, string D)
{
    {
        int x = 0;
        int y = 0;

        //------------------------------------------------------
        // Follow the plotmap line by line and make a sparse array in X and Y
        // of all the lines
        // 
        // Don't try to merge lines yet
        //------------------------------------------------------

        long tLoadArrays1 = GetTickCount();

        for (int i = 0; i < N; i++)
        {
            int li = L[i];
            char dc = D[i];

            // this maps the to direction.
            switch (dc)
            {
            case 'U':
            {
                pair<StartStopSet::iterator, bool> p = VStrokesRaw[x].emplace(StartStop(y, y + li));
                if (!p.second)
                {
                    p.first->stop = max(p.first->stop, y);
                }
                y += li;
                break;
            }
            case 'D':
            {
                pair<StartStopSet::iterator, bool> p = VStrokesRaw[x].emplace(StartStop(y - li, y));
                if (!p.second)
                {
                    p.first->stop = max(p.first->stop, y);
                }
                y -= li;
                break;
            }
            case 'L':
            {
                pair<YandXStopSet::iterator, bool> p = HStrokes[x-li].emplace(YandXStop(y, x)); // at X, there is a set, insert by Y
                if (!p.second)
                {
                    // already an entry at Y, extend the XStop
                    p.first->XStop = max(p.first->XStop, x);
                }
                x -= li;
                break;
            }
            case 'R':
            {
                pair<YandXStopSet::iterator, bool> p = HStrokes[x].emplace(YandXStop(y, x + li));
                if (!p.second)
                {
                    p.first->XStop = max(p.first->XStop, x+li);
                }
                x += li;
                break;
            }
            } // end switch
            if (x < xmin)
            {
                xmin = x;
            }
            if (x > xmax)
            {
                xmax = x;
            }
        }

        long tLoadArrays2 = GetTickCount() - tLoadArrays1;
        std::cout << "time for " << N << "entries, load arrays: " << tLoadArrays2 << endl;
    }

    if (VStrokesRaw.size() == 0 || HStrokes.size() == 0)
    {
        std::cout << "no lines." << endl;
        return 0;
    }

    long tGatherStrokes1 = GetTickCount();

    //-----------------------------------------------------------------
    // consolidate VStrokes
    //-----------------------------------------------------------------

    VStrokes = vector<X_YArray>(VStrokesRaw.size()); // have to make a copy since stuff in a set isn't changable normally.

    long tConsolidateArray1 = GetTickCount();
    long mergedVStrokes = 0;

    int xindex = 0;
    for (const pair<int, StartStopSet>& vs_kvp : VStrokesRaw)
    {
        int x = vs_kvp.first;
        const StartStopSet& vss = vs_kvp.second;
        int StrokeCount = vss.size();

        // set the initial size in the copy to be the same, but it will end up shorter if we merge any strokes
        VStrokes[xindex].Strokes.resize(StrokeCount);
        VStrokes[xindex].X = x;

        // go through each of the strokes at this X....
        int outStrokes = 0;
        StartStopSet::iterator i = vss.begin();
        int start = i->start;
        int stop = i->stop;
        VStrokes[xindex].Strokes[outStrokes].start = start;
        VStrokes[xindex].Strokes[outStrokes].stop = stop;
        outStrokes++;
        i++;
        for (; i != vss.end(); i++)
        {
            start = i->start;
            stop = i->stop;

            // and see if any of the existing merged strokes engulf this stroke and we can extend this one or throw it out.
            // we could probably use a binary search for this...
            int foundIndex = binarySearchMap2(VStrokes[xindex].Strokes, 0, VStrokes[xindex].Strokes.size() - 1, start);
            if (foundIndex == -1)
            {
                assert(outStrokes < VStrokes[xindex].Strokes.size());
                VStrokes[xindex].Strokes[outStrokes].start = start;
                VStrokes[xindex].Strokes[outStrokes].stop = stop;
                outStrokes++;
            }
            else
            {
                VStrokes[xindex].Strokes[foundIndex].stop = max(VStrokes[xindex].Strokes[foundIndex].stop, stop);
                mergedVStrokes++;
            }
        }
        VStrokes[xindex].Strokes.resize(outStrokes);

        xindex++;
    }
    VStrokesRaw.clear();

    long tConsolidateArray2 = GetTickCount() - tConsolidateArray1;;
    std::cout << "time to consolidate VStrokes " << tConsolidateArray2 << ", merged " << mergedVStrokes << " strokes" << endl;

    //-----------------------------------------------------------------
    // Run SweepX from all the way left to all the way right. We will step
    // to either the next new vertical stroke, or we'll step to the next place
    // an HStroke starts or stops. When an hstroke starts or stops, we'll update
    // a list (map) of "live lines" that represent the live horz strokes that are present
    // at SweepX. When SweepX comes to a new vertical stroke, the intersections are checked
    // one unit AFTER the vertical stroke's X location.
    //-----------------------------------------------------------------

    map<int, YandXStopSet>::iterator xiter = HStrokes.begin();
    map<int, YandXStopSet>::iterator xnext = xiter;
    xnext++;
    vector<X_YArray>::iterator vstrokeiter = VStrokes.begin();
    int vstrokex = vstrokeiter->X;

    long tFindPlusses1 = GetTickCount();

    int startX = xiter->first;
    int sweepX = startX;
    int checkRemoveCountdown = 0;
    int advancedX = 0;
    int xnextx = xnext->first;

    AddNewLinesToHLive(xiter);

    while(true)
    {
        if (xnext != HStrokes.end())
        {
            // only add new lines to Live list when we are about to cross a boundary
            // where there is a new stroke

            if (sweepX == xnextx - 1)
            {
                AddNewLinesToHLive(xnext);
            }

            if (sweepX == xnextx)
            {
                // time to swap iterators
                xiter = xnext;
                if (xnext != HStrokes.end())
                {
                    xnext++;
                    if (xnext != HStrokes.end())
                    {
                        xnextx = xnext->first;
                    }
                    else
                    {
                        xnextx = xmax;
                    }
                }
                checkRemoveCountdown = 2;
                advancedX++;
                if (advancedX % 1000 == 0)
                {
                    // cout << "advanced " << advancedX << " out of " << HStrokes.size() << ", hlive count = " << HLive.size() << endl;
                }
            }
        }

        // see if there are any possible vertical lines at this sweepX. if not, continue to
        // increment the iterator for the vertical strokes. Neat-o.

        while (vstrokex < sweepX )
        {
            vstrokeiter++;
            if (vstrokeiter == VStrokes.end())
            {
                break; // all done, no more vstrokes.
            }
            vstrokex = vstrokeiter->X;
            // do not advance xi or x to vstrokex, because it will mess up the 
            // optimized set of live horizontal lines at position vstrokex.
            // Just setting vstrokex ahead will optimize us to not try and find
            // plusses until we get there
        }

        if (vstrokeiter == VStrokes.end())
        {
            break;
        }

        if (vstrokex == sweepX)
        {
            // CheckForPlusses(vstrokeiter);5
        }

        //---------------------------------------------------------------------
        // remove old lines in HLive based on their xstop being less than current sweepX
        // Only remove lines one unit AFTER we've seen a new horz line
        //---------------------------------------------------------------------

        if (checkRemoveCountdown > 0)
        {
            checkRemoveCountdown--;
            if (checkRemoveCountdown == 0)
            {
                RemoveOldHLiveStrokes(sweepX);
            }
        }

        if (xnext == HStrokes.end() && HLive.size() == 0)
        {
            break;
        }

        sweepX++;
    }

    long tFindPlusses2 = GetTickCount() - tFindPlusses1;
    std::cout << "time for " << N << " entries, run sweep: " << tFindPlusses2 << endl;

    return foundPlusses.size();
}

void TryNPlusses(int N)
{
    int Li = 10000;
    vector<int> v(N);
    string s;
    for (int i = 0; i < N; i++)
    {
        v[i] = rand() % Li + 1;
        int dir = rand() % 4;
        if (i == 0) dir = 0;
        if (dir == 0) s += 'U';
       if (dir == 1) s += 'D';
        if (dir == 2) s += 'L';
        if (dir == 3) s += 'R';
    }
    long t1 = GetTickCount();
    long long plusses = getPlusSignCount(N, v, s);
    long t2 = GetTickCount() - t1;

    std::cout << "plusses for " << N << " entries = " << plusses << ", and took time " << t2 << endl;
}

void permute(string& sz, int i)
{
    // swap 1st with each character in the remainder

    for (int j = i ; j < sz.length(); j++)
    {
        char t = sz[i];
        sz[i] = sz[j];
        sz[j] = t;
        if (j < sz.length() - 1)
        {
            permute(sz, i + 1);
        }
        cout << sz << endl;
    }
}

int main()
{
    srand(1);
    int plusses;

    #if false
    plusses = getPlusSignCount(8, { 1, 2, 2, 1, 1, 2, 2, 1 }, "UDUDLRLR");
    std::cout << "ex 3: " << plusses << endl; // should be 1
    #endif

#if false
    plusses = getPlusSignCount(9, { 6, 3, 4, 5, 1, 6, 3, 3, 4 }, "ULDRULURD");
    std::cout << "ex 1: " << plusses << endl; // should be 4
#endif

#if false
    plusses = getPlusSignCount(8, { 1, 1, 1, 1, 1, 1, 1, 1 }, "RDLUULDR");
    std::cout << "ex 2: " << plusses << endl;
#endif

    // TryNPlusses(20);
    // TryNPlusses(200);
    // TryNPlusses(2000);
    // TryNPlusses(20000);
    // TryNPlusses(200000);
    TryNPlusses(2000000);
    TryNPlusses(4000000);
    TryNPlusses(8000000);
    cout << maxLinesInHLive << endl;

}
