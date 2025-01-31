
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PBRT_CORE_KDTREE_H
#define PBRT_CORE_KDTREE_H

// core/kdtree.h*
#include "pbrt.h"
#include "geometry.h"

// KdTree Declarations
struct KdNode {
    void init(const mfloat<length_d> &p, uint32_t a) {
        splitPos = p;
        splitAxis = a;
        rightChild = (1<<29)-1;
        hasLeftChild = 0;
    }
    void initLeaf() {
        splitAxis = 3;
        rightChild = (1<<29)-1;
        hasLeftChild = 0;
    }
    // KdNode Data
    mfloat<length_d> splitPos;
    uint32_t splitAxis:2;
    uint32_t hasLeftChild:1;
    uint32_t rightChild:29;
};


template <typename NodeData> class KdTree {
public:
    // KdTree Public Methods
    KdTree(const std::vector<NodeData> &data);
    ~KdTree() {
        FreeAligned(nodes);
        delete[] nodeData;
    }
    template <typename LookupProc, typename S> void Lookup(const Point<S> &p,
            LookupProc &process, mfloat<area_d> &maxDistSquared) const;
private:
    // KdTree Private Methods
    void recursiveBuild(uint32_t nodeNum, int start, int end,
        const NodeData **buildNodes);
    template <typename LookupProc, typename S> void privateLookup(uint32_t nodeNum,
        const Point<S> &p, LookupProc &process, mfloat<area_d> &maxDistSquared) const;

    // KdTree Private Data
    KdNode *nodes;
    NodeData *nodeData;
    uint32_t nNodes, nextFreeNode;
};


template <typename NodeData> struct CompareNode {
    CompareNode(int a) { axis = a; }
    int axis;
    bool operator()(const NodeData *d1, const NodeData *d2) const {
        return d1->p[axis] == d2->p[axis] ? (d1 < d2) :
                                            d1->p[axis] < d2->p[axis];
    }
};



// KdTree Method Definitions
template <typename NodeData>
KdTree<NodeData>::KdTree(const std::vector<NodeData> &d) {
    nNodes = d.size();
    nextFreeNode = 1;
    nodes = AllocAligned<KdNode>(nNodes);
    nodeData = new NodeData[nNodes];
    std::vector<const NodeData *> buildNodes(nNodes, NULL);
    for (uint32_t i = 0; i < nNodes; ++i)
        buildNodes[i] = &d[i];
    // Begin the KdTree building process
    recursiveBuild(0, 0, nNodes, &buildNodes[0]);
}


template <typename NodeData> void
KdTree<NodeData>::recursiveBuild(uint32_t nodeNum, int start, int end,
        const NodeData **buildNodes) {
    // Create leaf node of kd-tree if we've reached the bottom
    if (start + 1 == end) {
        nodes[nodeNum].initLeaf();
        nodeData[nodeNum] = *buildNodes[start];
        return;
    }

    // Choose split direction and partition data

    // Compute bounds of data from _start_ to _end_
    BBox<world_s> bound;
    for (int i = start; i < end; ++i)
        bound = Union(bound, buildNodes[i]->p);
    int splitAxis = bound.MaximumExtent();
    int splitPos = (start+end)/2;
    std::nth_element(&buildNodes[start], &buildNodes[splitPos],
                     &buildNodes[end], CompareNode<NodeData>(splitAxis));

    // Allocate kd-tree node and continue recursively
    nodes[nodeNum].init(buildNodes[splitPos]->p[splitAxis], splitAxis);
    nodeData[nodeNum] = *buildNodes[splitPos];
    if (start < splitPos) {
        nodes[nodeNum].hasLeftChild = 1;
        uint32_t childNum = nextFreeNode++;
        recursiveBuild(childNum, start, splitPos, buildNodes);
    }
    if (splitPos+1 < end) {
        nodes[nodeNum].rightChild = nextFreeNode++;
        recursiveBuild(nodes[nodeNum].rightChild, splitPos+1,
                       end, buildNodes);
    }
}


template <typename NodeData> template <typename LookupProc, typename S>
void KdTree<NodeData>::Lookup(const Point<S> &p, LookupProc &proc,
        mfloat<area_d> &maxDistSquared) const {
    privateLookup(0, p, proc, maxDistSquared);
}


template <typename NodeData> template <typename LookupProc, typename S>
void KdTree<NodeData>::privateLookup(uint32_t nodeNum, const Point<S> &p,
        LookupProc &process, mfloat<area_d> &maxDistSquared) const {
    KdNode *node = &nodes[nodeNum];
    // Process kd-tree node's children
    int axis = node->splitAxis;
    if (axis != 3) {
        mfloat<area_d> dist2 = (p[axis] - node->splitPos) * (p[axis] - node->splitPos);
        if (p[axis] <= node->splitPos) {
            if (node->hasLeftChild)
                privateLookup(nodeNum+1, p, process, maxDistSquared);
            if (dist2 < maxDistSquared && node->rightChild < nNodes)
                privateLookup(node->rightChild, p, process, maxDistSquared);
        }
        else {
            if (node->rightChild < nNodes)
                privateLookup(node->rightChild, p, process, maxDistSquared);
            if (dist2 < maxDistSquared && node->hasLeftChild)
                privateLookup(nodeNum+1, p, process, maxDistSquared);
        }
    }

    // Hand kd-tree node to processing function
    mfloat<area_d> dist2 = DistanceSquared(nodeData[nodeNum].p, p);
    if (dist2 < maxDistSquared)
        process(p, nodeData[nodeNum], dist2, maxDistSquared);
}



#endif // PBRT_CORE_KDTREE_H
