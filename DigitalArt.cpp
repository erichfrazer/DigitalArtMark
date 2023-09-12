#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <windows.h>
#include <assert.h>
#include <string>
using namespace std;

int tries = 0;

// Point
struct Point
{
    int x;
    int y;
};

// Segment
struct Segment
{
    // saves memory by only having an index and start/stop.
    // for vertical lines, the index is Y, and start/end refer to X range
    // for horizontal lines, the index is X, and start/end refer to Y range
    int index;
    int start;
    int end;

    Segment( int i, int s, int e )
    {
        index = i;
        start = s;
        end = e;
    }

    Segment()
    {
        // needed for resize()! Dont' do anything
    }
};

// for hLive list of horz strokes
struct StartStop
{
    int start;
    int stop;

    StartStop(int s, int e)
    {
        start = s;
        stop = e;
    }

    StartStop()
    {
        start = 0;
        stop = 0;
    }
};

class CrossDetector
{
    vector<Segment> segHorizontal;
    vector<Segment> segVertical;
    map<int, StartStop> hLive;

    void AddMovesToPool(int N, vector<int> L, string D);
    void SortPoolForDecimation();
    void DecimatePool();
    void SortPoolForCounting();
    int CountCrosses();

    static int SortHorizontal(const void* item1, const void* item2);
    static int SortVertical(const void* item1, const void* item2);
    static int SortHorizontalToVertical(const void* item1, const void* item2);

public:

    CrossDetector()
    {
        segHorizontal.clear();
        segVertical.clear();
        hLive.clear();
    }

    long long getPlusSignCount(int N, vector<int> L, string D);

};

// Simple move
struct Move
{
    char direction;  // N S E W
    int distance;    // Segment length
};

// Add moves to segment pool
void CrossDetector::AddMovesToPool(int N, vector<int> L, string D)
{
    tries++;

    // We start at the origin
    Point ptCurrent = { 0, 0 };

    for (int i = 0; i < N; ++i)
    {
        // Each move creates a segment that starts where we are, and ends where we end
        Point ptEnd = ptCurrent;

        // Add the segment to the pool
        // Ensure we are always top to bottom, left to right for optimization later
        int distance = L[i];
        switch (D[i])
        {
        case 'U':
            ptEnd.y += distance;
            segVertical.push_back(Segment(ptCurrent.x, ptCurrent.y, ptEnd.y));
            break;
        case 'D':
            ptEnd.y -= distance;
            segVertical.push_back(Segment(ptCurrent.x, ptEnd.y, ptCurrent.y));
            break;
        case 'R':
            ptEnd.x += distance;
            segHorizontal.push_back(Segment(ptCurrent.y, ptCurrent.x, ptEnd.x));
            break;
        case 'L':
            ptEnd.x -= distance;
            segHorizontal.push_back(Segment(ptCurrent.y, ptEnd.x, ptCurrent.x));
            break;
        }

        // End point is now current point
        ptCurrent = ptEnd;
    }
}

// For horizontal segments, sort by lowest Y, then lowest X, then length (small to large)
int CrossDetector::SortHorizontal(const void* item1, const void* item2)
{
    const Segment* seg1 = (Segment*)item1;
    const Segment* seg2 = (Segment*)item2;

    // Sort on Y
    if (seg1->index < seg2->index)
    {
        return -1; // seg2 is more on Y
    }
    else if (seg1->index > seg2->index)
    {
        return 1;  // seg1 is more on Y
    }
    else
    {
        // Next level sort (X)
        if (seg1->start < seg2->start)
        {
            return -1;  // seg1 is lower on X
        }
        else if (seg1->start > seg2->start)
        {
            return 1;   // seg2 is lower on X
        }
        else
        {
            // Next level sort (length)
            if (seg1->end < seg2->end)
            {
                return -1; // seg1 is shorter
            }
            else if (seg1->end > seg2->end)
            {
                return 1;  // seg2 is shorter
            }
            else
            {
                return 0; // Same segment
            }
        }
    }

    // Should not get here
    return 0;
}


// For vertical segments, sort by smallest X, then lowest Y, then length (small to large)
int CrossDetector::SortVertical(const void* item1, const void* item2)
{
    const Segment* seg1 = (Segment*)item1;
    const Segment* seg2 = (Segment*)item2;

    // Sort on X
    if (seg1->index < seg2->index)
    {
        return -1; // seg1 is less on X
    }
    else if (seg1->index > seg2->index)
    {
        return 1;  // seg2 is less on X
    }
    else
    {
        // Next level sort (Y)
        if (seg1->start < seg2->start)
        {
            return -1;  // seg2 is higher on Y
        }
        else if (seg1->start > seg2->start)
        {
            return 1;   // seg1 is higher on Y
        }
        else
        {
            // Next level sort (length)
            if (seg1->end > seg2->end)
            {
                return -1; // seg1 is shorter
            }
            else if (seg1->end < seg2->end)
            {
                return 1;  // seg2 is shorter
            }
            else
            {
                return 0; // Same segment
            }
        }
    }

    // Should not get here
    return 0;
}

// resort the horizontal lines to be ordered by increasing X, then increasing Y
int CrossDetector::SortHorizontalToVertical(const void* item1, const void* item2)
{
    const Segment* seg1 = (Segment*)item1;
    const Segment* seg2 = (Segment*)item2;

    // Sort on X
    if (seg1->start < seg2->start)
    {
        return -1; // seg1 is less on X
    }
    else if (seg1->start > seg2->start)
    {
        return 1;  // seg2 is less on X
    }
    else
    {
        // Next level sort (Y)
        if (seg1->index < seg2->index)
        {
            return -1;  // seg2 is higher on Y
        }
        else if (seg1->index > seg2->index)
        {
            return 1;   // seg1 is higher on Y
        }
        else
        {
            // Next level sort (length)
            if (seg1->end > seg2->end)
            {
                return -1; // seg1 is shorter
            }
            else if (seg1->end < seg2->end)
            {
                return 1;  // seg2 is shorter
            }
            else
            {
                return 0; // Same segment
            }
        }
    }

    // Should not get here
    return 0;
}

// Sort the pool
// Segments already in x min to x max, and y max to y min
void CrossDetector::SortPoolForDecimation()
{
    qsort(segHorizontal.data(), segHorizontal.size(), sizeof(Segment), SortHorizontal);
    qsort(segVertical.data(), segVertical.size(), sizeof(Segment), SortVertical);
}

void CrossDetector::DecimatePool()
{
    // Loop through and combine segments with the following properties
    // Given they are sorted, we can likely combine them

    if (tries == 5)
    {
        int stop = 0;
    }

    // First vertical segements
    {
        int iFinalVertical = 0;
        for (int iSegVertical = 1; iSegVertical < segVertical.size(); ++iSegVertical)
        {
            // Check if current segment can be collapse into final set

            // We assume sorted by x, top to bottom vertical segements

            // If not on the same X, nothing to do
            if (segVertical[iFinalVertical].index != segVertical[iSegVertical].index)
            {
                // Commit this one to the final list
                if (segVertical[iFinalVertical].end - segVertical[iFinalVertical].start > 1)
                {
                    ++iFinalVertical;
                }
                segVertical[iFinalVertical] = segVertical[iSegVertical];
            }
            else
            {
                // On the same X, see if segment overlaps or extends
                if (segVertical[iSegVertical].start <= segVertical[iFinalVertical].end)
                {
                    // Overlap, take the furthest end
                    if (segVertical[iSegVertical].end > segVertical[iFinalVertical].end)
                    {
                        segVertical[iFinalVertical].end = segVertical[iSegVertical].end;
                    }
                }
                else
                {
                    // Done - Check to ensure we have a larger than one length segement (otherwise it can't create a cross)
                    if (segVertical[iFinalVertical].end - segVertical[iFinalVertical].start > 1)
                    {
                        ++iFinalVertical;
                    }
                    segVertical[iFinalVertical] = segVertical[iSegVertical];
                }
            }
        }

        // Only take last segment if greater than one in length
        if (segVertical[iFinalVertical].end - segVertical[iFinalVertical].start > 1)
        {
            ++iFinalVertical;
        }
        int removed = segVertical.size() - iFinalVertical;
        // cout << "removing " << removed << " vertical lines, lines remaining: " << iFinalVertical << endl;
        segVertical.resize(iFinalVertical);
    }

    {
        // Second horizontal segements
        int iFinalHorizontal = 0;
        for (int iSegHorizontal = 0; iSegHorizontal < segHorizontal.size(); ++iSegHorizontal)
        {
            // Check if current segment can be collapsed into final set

            // We assume sorted by y, then by length

            // If not on the same y, nothing to do
            if (segHorizontal[iFinalHorizontal].index != segHorizontal[iSegHorizontal].index)
            {
                // Commit this one to final list if it is longer than length one
                if (segHorizontal[iFinalHorizontal].end - segHorizontal[iFinalHorizontal].start > 1)
                {
                    ++iFinalHorizontal;
                }
                segHorizontal[iFinalHorizontal] = segHorizontal[iSegHorizontal];
            }
            else
            {
                // On the same Y, see if segment overlaps or extends
                if (segHorizontal[iSegHorizontal].start <= segHorizontal[iFinalHorizontal].end)
                {
                    // Overlap, take the furthest end
                    if (segHorizontal[iSegHorizontal].end > segHorizontal[iFinalHorizontal].end)
                    {
                        segHorizontal[iFinalHorizontal].end = segHorizontal[iSegHorizontal].end;
                    }
                }
                else
                {
                    // Done - Check to ensure we have a larger than one length segement (otherwise it can't create a cross)
                    if (segHorizontal[iFinalHorizontal].end - segHorizontal[iFinalHorizontal].start > 1)
                    {
                        ++iFinalHorizontal;
                    }
                    segHorizontal[iFinalHorizontal] = segHorizontal[iSegHorizontal];
                }
            }

        }

        // Only take last segment if greater than one in length
        if (segHorizontal[iFinalHorizontal].end - segHorizontal[iFinalHorizontal].start > 1)
        {
            ++iFinalHorizontal;
        }
        int removed = segHorizontal.size() - iFinalHorizontal;
        // cout << "removing " << removed << " horizontal lines, lines remaining: " << iFinalHorizontal << endl;
        segHorizontal.resize(iFinalHorizontal);
    }
}

// Sort pool on X, then Y, then length
void CrossDetector::SortPoolForCounting()
{
    // Vertical is already sorted by min X to max X, min Y to max Y, then by length
    // Sort horizontal similar to vertical to speed up final cross check
    // original segHorizontal: index is by Y, and start and stop are by X.
    qsort(segHorizontal.data(), segHorizontal.size(), sizeof(Segment), SortVertical);
}

int CrossDetector::CountCrosses()
{
    int nCrosses = 0;

    // We assume segments come in sorted by X min to X max, then Y max to Y min, then by length
    // As such, we can scan through only checking those that likely will intersect (effectively moving diagonally)

    // pretty sure this is buggy
#if false
    // For each vertical, only check the horizontals that could intersect
    int iHorizontalStart = 0;
    for (int iVertical = 0; iVertical < nVertical; ++iVertical)
    {`
        for (int iHorizontal = iHorizontalStart; iHorizontal < nHorizontal; ++iHorizontal)
        {
            // Check for segment overlapp
            if (segHorizontal[iHorizontal].start.x >= segVertical[iVertical].end.x)
            {
                // Past all horizontal segments that matter, next vertical
                break;
            }

            if (segHorizontal[iHorizontal].end.x <= segVertical[iVertical].end.x)
            {
                // Past all horizontal segements that matter, in this case, don't look at those again
                ++iHorizontalStart;
            }
            else if (segHorizontal[iHorizontal].start.y < segVertical[iVertical].start.y &&
                segHorizontal[iHorizontal].end.y > segVertical[iVertical].end.y)
            {
                ++nCrosses;
            }
        }
    }
#endif

    /*
    *   <--    ----    ----    ----    ----    ----    ---
    * >---|----|--|----|--|----|--|----|--|----|--|----|--|---
    *     ------  ------  ------  ------  ------  ------  ---|
    * problem 1: could be a line from very beginning that extends all the way
    * to the right.
    *
    * ----------------------------|---------------------------
    *                             |                          |
    * ----------------------------|---------------------------
    * |                           |
    * ----------------------------|---------------------------
    *                             |                          |
    * ----------------------------|---------------------------
    * |                           |
    * -----------------------------
    * problem 2: bunch of horizontal lines then a vertical one stabbing thru them
    *
    * problem 3: combination of 1 and 2. visualize it!
    */

    int vix = 0;
    int hix = 0;
    
    map<int, StartStop>::iterator hiter = hLive.begin();

    // caching variables so we don't have to do lookups, and we test against these a lot from loop to loop
    int vx = segVertical[vix].index;
    int vsy = segVertical[vix].start;
    int vey = segVertical[vix].end;
    int vix_limit = segVertical.size();
    int hix_limit = segHorizontal.size();
    bool bFirstLine = true;

#define hStrokeY(x) x->first

     while (true)
     {
         if (vix == vix_limit)
         {
             // cout << "ran out of vertical strokes, we're done" << endl;
             break;
         }

        
         if (bFirstLine)
         {
            bFirstLine = false;
            goto continue_y;
         }

        //----------------------------------------------
        // Advance vix as long as vstroke is out of range (less than) of hiter's Y.
        // It might advance to the next vstroke, which is fine
        //----------------------------------------------

        if(hiter == hLive.end() || vey <= hStrokeY(hiter))
        {
            vix++;
            if (vix == vix_limit)
            {
                // cout << "Done, because we needed a new vstroke, but were at the end of vstrokes" << endl;
                break;
            }

            // cout << "viter++, because ran out of hLines, or our stop endpoint " << vey << " is less than hstrokeY " << endl;
            vx = segVertical[vix].index;
            vsy = segVertical[vix].start;
            vey = segVertical[vix].end;
            goto continue_y;
        }

        //----------------------------------------------
        // Advance hLive iter as long as horz stroke's Y is <= the current stroke's start
        // It might go off the end of the list, which is okay.
        //----------------------------------------------

        if (hStrokeY(hiter) <= vsy)
        {
            // cout << "hiter++, because hStrokeY @ " << hStrokeY(hiter) << " <= vstroke's StartY @ " << vsy << endl;
            hiter++;
            continue; // not continue_y
        }

        //----------------------------------------------
        // hiter's Y > vstroke's StartY, so it's a candidate for a plus.
        // But vstroke's StopY could also be less than hStrokeY.
        // Check for a plus at this hStrokeY
        //----------------------------------------------

        if (hStrokeY(hiter) < vey && hiter->second.start < vx)
        {
            // cout << "* Found a cross at " << vx << "," << hStrokeY(hiter) << endl;
            nCrosses++;
        }

        // cout << "hiter++" << endl;

        hiter++;
        continue; // not continue_y

    continue_y:

        //----------------------------------------------
        // Add horz lines to hLive. Only do this after checking
        // for plusses.
        //----------------------------------------------

        if (vix < vix_limit)
        {
            // if we're on the last vstroke on this vertical line and have already
            // checked for plusses

            if (vix == 0 || vx != segVertical[vix-1].index)
            {
                // cout << "Moving to a new vstroke " << vix << " at VX = " << vx << ", vstroke goes from Y=" << segVertical[vix].start.y << " to " << segVertical[vix].end.y << endl;

                //----------------------------------------------
                // If at new vstroke, remove all dead hLive lines
                //----------------------------------------------

                hiter = hLive.begin();
                while (hiter != hLive.end())
                {
                    if (hiter->second.stop <= vx)
                    {
                        // cout << "  at VX=" << vx << " removing hLive line @ Y = " << hiter->first << " from " << hiter->second.start << " to " << hiter->second.stop << endl;
                        hiter = hLive.erase(hiter);
                    }
                    else
                    {
                        hiter++;
                    }
                }

                //----------------------------------------------
                // Add new hstrokes to hLive
                //----------------------------------------------

                if (hix != hix_limit)
                {
                    // after sorting horz by X, the index is
                    while (segHorizontal[hix].start <= vx)
                    {
                        int startX = segHorizontal[hix].start;
                        int stopX = segHorizontal[hix].end;
                        int horzY = segHorizontal[hix].index;
                        hLive[horzY] = StartStop(startX, stopX);
                        // cout << "  at VX=" << vx << " adding hLive line @ Y = " << horzY << " from " << startX << " to " << stopX << endl;
                        hix++;
                        if (hix == hix_limit)
                        {
                            break;
                        }
                    }
                }

                if (hLive.size() != 0)
                {
                    hiter = hLive.begin();
                }
                else
                {
                    // cout << "Ran out of hLive lines, we must be done?" << endl;
                    break; // no lines, we must be done
                }
            }
        }
        else
        {
            assert(false);
        }

    } // every vertical stroke

    return nCrosses;
}

// Calculate cross count
long long CrossDetector::getPlusSignCount(int N, vector<int> L, string D) 
{
    int nCrosses = 0;

    // Add moves to segment pool
    AddMovesToPool(N, L, D);

    // Sort the pool
    SortPoolForDecimation();

    // Decimate pool
    DecimatePool();

    // Sort the pool for counting
    SortPoolForCounting();

    // Count Crosses
    nCrosses = CountCrosses();

    return nCrosses;
}

long long testN(int n, vector<int> v = vector<int>(), string dir = string())
{
    long t1 = GetTickCount();
    for (int i = 0; i < n; i++)
    {
        int len = 1 + (rand() % 200);
        v.push_back(len);
        int c = rand() % 4;
        if (c == 0) dir += "U";
        if (c == 1) dir += "D";
        if (c == 2) dir += "L";
        if (c == 3) dir += "R";
    }
    CrossDetector cd;
    long long result = cd.getPlusSignCount(n, v, dir);
    long t2 = GetTickCount();
    printf("Time for %d entries: %d\n\n\n\n", n, t2 - t1);
    return result;
}

int main()
{
    // vector<int> v = { 1, 2, 2, 1, 1, 2, 2, 1 };
    // string dir = "UDUDLRLR";

#if false
    vector<int> v = { 6, 3, 4, 5, 1, 6, 3, 3, 4 };
    string dir = "ULDRULURD";
    long long result = testN(9, v, dir);
#endif

    for (int i = 0 ; i < 1 ; i++ )
    {
        long long ll = testN(200);
    }

    return -1;

    for (int i = 32; i < 200000; i++)
    {
        long long ll = testN(i);
    }
}

