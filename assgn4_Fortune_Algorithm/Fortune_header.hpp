// NAME:Akash Tiwari
// ROLL:17CS10003

#include "point.hpp"
#include "DCEL.hpp"

#ifndef Voronoi_Structure_hpp
#define Voronoi_Structure_hpp

#define POINT_EPSILON 1.0e-6

using namespace std;

// Beachline headers
class Event;
namespace beachline 
    {   
    using namespace DCEL;
    class BL_Node;
    typedef shared_ptr<BL_Node> BL_NPtr;
    class BL_Node 
    {
    public:
        int height;
        double *sweepline;
        const vector<Pointxy> *points;
        pair<int, int> indices;
        BL_NPtr left, right, parent;
        shared_ptr<Event> circle_event;
        shared_ptr<H_edge> edge;
        BL_Node(const pair<int,int>& _indices, double* _sweepline = nullptr, const vector<Pointxy>* _points = nullptr, BL_NPtr _left = nullptr, BL_NPtr _right = nullptr, BL_NPtr _parent = nullptr, int _height = 1);
        BL_NPtr next, prev;    
        inline bool is_leaf() 
        {
            return indices.first == indices.second;
        }
        inline int get_id() 
        {
            return indices.first;
        }
        inline bool has_indices(int a, int b) 
        {
            return indices.first == a && indices.second == b;
        }
        inline bool has_indices(const pair<int,int> &p) 
        {
            return indices.first == p.first && indices.second == p.second;
        }
        double value();        
    };
    void connect(BL_NPtr prev, BL_NPtr next);
    bool is_root(BL_NPtr node);
    // find the height of a node
    int find_h(BL_NPtr node);
    // change the height of a node
    void change_h(BL_NPtr node);
    // rotate left around the node
    BL_NPtr rotate_left(BL_NPtr node);
    // rotate right around the node
    BL_NPtr rotate_right(BL_NPtr node);
    // finds a leaf node in tree such that is under the parabolic arc corresponding to this leaf
    BL_NPtr find(BL_NPtr root, double x);
    // replace a leaf node with a new subtree
    BL_NPtr replace(BL_NPtr node, BL_NPtr new_node);
    // find balance of the node
    int fetch_balance(BL_NPtr node);
    BL_NPtr remove(BL_NPtr leaf);;
    pair<BL_NPtr, BL_NPtr> breakpoints(BL_NPtr leaf);
    BL_NPtr make_subtree(int index, int index_behind, double *sweepline, const vector<Pointxy> *points, vector<H_EdgePtr> &edges);
    BL_NPtr make_simple_subtree(int index, int index_behind, double *sweepline, const vector<Pointxy> *points, vector<H_EdgePtr> &edges);  
    bool _validate(BL_NPtr node);
    bool _check_balance(BL_NPtr node);
    void print_tree(BL_NPtr root, int width = 7);   
}
namespace bl = beachline;
void build_voronoi(const vector<Pointxy> &points, vector<bl::H_EdgePtr> &halfedges, vector<bl::VertexPtr> &vertices, vector<bl::H_EdgePtr> &faces);
#endif
