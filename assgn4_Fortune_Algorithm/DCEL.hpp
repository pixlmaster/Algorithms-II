// NAME:Akash Tiwari
// ROLL:17CS10003


#include "point.hpp"

#ifndef DCEL_hpp
#define DCEl_hpp

#define POINT_EPSILON 1.0e-6

using namespace std;

namespace DCEL 
{
    class Vertex;
    class H_edge;
    typedef shared_ptr<H_edge> H_EdgePtr;
    typedef shared_ptr<Vertex> VertexPtr;
    class Vertex 
    {
    public:   
        Pointxy point;
        H_EdgePtr edge;
        Vertex(const Pointxy &pos, H_EdgePtr incident_edge = nullptr);
        inline double x()
        {
            return point.x; 
        }
        inline double y() 
        {
            return point.y; 
        }
    };

    class H_edge 
    {
    public:   
        int l_index, r_index;
        VertexPtr vertex;
        H_EdgePtr twin;
        H_EdgePtr next;
        H_EdgePtr prev;
        H_edge(int _l_index, int _r_index, VertexPtr _vertex = nullptr);
        inline VertexPtr vertex0()
        {
            return vertex; 
        }
        inline VertexPtr vertex1()
        {
            return twin->vertex; 
        }
        inline bool is_finite() 
        {
            return vertex != nullptr && twin->vertex != nullptr;
        }
        // Iterators around vertex
        H_EdgePtr vertexNextCCW();
        H_EdgePtr vertexNextCW();
    };

    pair<H_EdgePtr, H_EdgePtr> make_twins(int left_index, int right_index);
    pair<H_EdgePtr, H_EdgePtr> make_twins(const pair<int,int> &indices);    
    void connect_h_edges(H_EdgePtr p1, H_EdgePtr p2);
}

#endif