// NAME:Akash Tiwari
// ROLL:17CS10003

#include "Fortune_header.hpp"
#define _DEBUG_
#define BREAKPOINTS_EPSILON 1.0e-5
#define CIRCLE_CENTER_EPSILON 1.0e-4
#if defined(_WIN64) || defined(_WIN32)      
    #define isnan(x) _isnan(x)
#endif

using namespace std;

#include <cmath>
#include <queue>


// functions for parabola
int intersectionPoints(const Pointxy &point_f1, const Pointxy &point_f2, double directrix) 
{
    if (fabs(point_f1.x - point_f2.x) < POINT_EPSILON && fabs(point_f1.y - point_f2.y) < POINT_EPSILON) 
    {
        return -1;
    }
    else if (fabs(point_f1.y - point_f2.y) < POINT_EPSILON)
    {
        return 1;
    }
    else
    {
        return 2;
    }
}

// finds the list of intersection point
vector<Pointxy> findIntersectionPoints(const Pointxy &point_f1, const Pointxy &point_f2, double dist) 
{
    vector<Pointxy> answer;
    if (fabs(point_f1.x - point_f2.x) < POINT_EPSILON) 
    {
        double y = 0.5*(point_f1.y+point_f2.y),D=sqrt(dist*dist-dist*(point_f1.y+point_f2.y)+point_f1.y*point_f2.y);
        answer.push_back(Pointxy(point_f1.x-D,y));
        answer.push_back(Pointxy(point_f1.x+D,y));
    } 
    else if (fabs(point_f1.y-point_f2.y) < POINT_EPSILON) 
    { 
        double x = 0.5*(point_f1.x+point_f2.x);
        answer.push_back(Pointxy(x, 0.5*((x-point_f1.x)*(x-point_f1.x)+point_f1.y*point_f1.y-dist*dist)/(point_f1.y-dist)));
    } 
    else 
    { 
        double D= 2. * sqrt(pow(point_f1.x-point_f2.x, 2)*(dist-point_f1.y)*(dist-point_f2.y)*(pow(point_f1.x-point_f2.x, 2) + pow(point_f1.y-point_f2.y, 2)));
        double doub_T= -2. * dist*pow(point_f1.x-point_f2.x, 2)+(point_f1.y+point_f2.y)*(pow(point_f2.x-point_f1.x, 2)+pow(point_f2.y-point_f1.y, 2));
        double doub_Q= 2. * pow(point_f1.y-point_f2.y, 2);   
        double y1=(doub_T-D)/doub_Q, y2=(doub_T+D)/doub_Q;
        double x1 = 0.5*(point_f1.x*point_f1.x-point_f2.x*point_f2.x+(2*y1-point_f2.y-point_f1.y)*(point_f2.y-point_f1.y))/(point_f1.x-point_f2.x);
        double x2 = 0.5*(point_f1.x*point_f1.x-point_f2.x*point_f2.x+(2*y2-point_f2.y-point_f1.y)*(point_f2.y-point_f1.y))/(point_f1.x-point_f2.x);
        if (x1<=x2) {} else
        {
            swap(x1, x2);
            swap(y1, y2); 
        }
        answer.push_back(Pointxy(x1, y1));
        answer.push_back(Pointxy(x2, y2));
    }
    return answer;
}

// 2-d points func definitions

const double Pointxy::Inf = numeric_limits<double>::infinity();
Pointxy::Pointxy_Compare Pointxy::xy_compare = Pointxy::Pointxy_Compare();
Pointxy::Pointxy(double _x, double _y) : x(_x), y(_y) 
{}
Pointxy::Pointxy(const Pointxy &point) : x(point.x), y(point.y) ,id(point.id) 
{}

Pointxy operator+(const Pointxy &point1, const Pointxy &point2) 
{
    return Pointxy(point1.x + point2.x, point1.y + point2.y);
}


ostream &operator<<(ostream &stream, const Pointxy &p) 
{
    stream << "(" << p.x << "," << p.y << ")";
    return stream;
}
vector<Pointxy> &operator<<(vector<Pointxy> &v, const Pointxy &p) 
{
    v.push_back(p);
    return v;
}
bool Pointxy::checkLeftTurn(const Pointxy &point1, const Pointxy &point2, const Pointxy &point3) 
{
    return (vector_prod(point2 - point1, point3 - point2) > 0.0);
}
void Pointxy::setY(double y) 
{
    this->y = y;
}
bool Pointxy::checkRightTurn(const Pointxy &point1, const Pointxy &point2, const Pointxy &point3) 
{
    return (vector_prod(point2 - point1, point3 - point2) < 0.0);
}
Pointxy operator/(const Pointxy &p, double value) 
{
    return Pointxy(p.x / value, p.y / value);
}
Pointxy operator-(const Pointxy &p) 
{
    return Pointxy(-p.x, -p.y);
}
double scalar_prod(const Pointxy &point1, const Pointxy &point2) 
{
    return point1.x * point2.x + point1.y * point2.y;
}
double vector_prod(const Pointxy &point1, const Pointxy &point2) 
{
    return point1.x * point2.y - point1.y * point2.x;
}
bool equal(const Pointxy &point1, const Pointxy &point2, double EPSILON) 
{
    return (fabs(point1.x - point2.x) < EPSILON && fabs(point1.y - point2.y) < EPSILON);
}
Pointxy &Pointxy::operator+=(const Pointxy &p) 
{
    x += p.x;
    y += p.y;
    return *this;
}
bool equal(double v1, double v2, double EPSILON) 
{
    return fabs(v1 - v2) < EPSILON;
}
Pointxy &Pointxy::operator-=(const Pointxy &p) 
{
    x -= p.x;
    y -= p.y;
    return *this;
}
Pointxy &Pointxy::operator*=(double value) 
{
    x *= value;
    y *= value;
    return *this;
}
void Pointxy::normalize() 
{
    double n = norm();
    x /= n;
    y /= n;
}
Pointxy operator*(const Pointxy &p, double value) 
{
    return Pointxy(p.x * value, p.y * value);
}
Pointxy operator*(double value, const Pointxy &p) 
{
    return Pointxy(p.x * value, p.y * value);
}
double Pointxy::operator[](int i) 
{
    if (i!=0) 
    {
        return y;
    }
    else
    {
        return x;
    }
}
void Pointxy::setX(double x) 
{
    this->x = x;
}

Pointxy &Pointxy::operator/=(double value) 
{
    x /= value;
    y /= value;
    return *this;
}
void Pointxy::setId(int id)
{
    this->id = id;

}
Pointxy operator-(const Pointxy &point1, const Pointxy &point2) 
{
    return Pointxy(point1.x - point2.x, point1.y - point2.y);
}
Pointxy operator/(const Pointxy &point1, const Pointxy &point2) 
{
    return Pointxy(point1.x / point2.x, point1.y / point2.y);
}
    
bool Pointxy::checkVertical() 
{
    return (y == Inf && !isnan(x) && x != Inf);
}
bool Pointxy::checkHorizontal() 
{
    return (x == Inf && !isnan(y) && y != Inf);
}
bool Pointxy::isValid() 
{
    if (x != Inf || y != Inf)
    {
        return (!isnan(x) && !isnan(y));
    }
    else
    {
        return false;
    }
}
Pointxy Pointxy::normalized() 
{
    return (*this) / this->norm();
}

double Pointxy::norm() 
{
    return sqrt(x * x + y * y);
}
double Pointxy::norm2() 
{
    return x *x + y * y;
}
Pointxy Pointxy::Rotate90CW() 
{
    return Pointxy(y, -x);
}
Pointxy Pointxy::Rotate90CCW() 
{
    return Pointxy(-y, x);
}

// func for beachline
// handles addition, deletion of leaf and re-balancing
namespace bl = beachline;
namespace beachline 
{
   BL_Node::BL_Node(const pair<int,int>& _indices,double* _sweepline ,const vector<Pointxy>* _points,BL_NPtr _left, BL_NPtr _right,BL_NPtr _parent,int _height) :
                                 indices(_indices), left(_left), right(_right),
                                  parent(_parent), height(_height),
                                  sweepline(_sweepline), points(_points),
                                  next(nullptr), prev(nullptr) {}
   
    // connect as a list
    void connect(BL_NPtr prev, BL_NPtr next) 
    {
        prev->next = next;
        next->prev = prev;
    }
    // find whether node is root or not
    bool is_root(BL_NPtr BL_node) 
    {
        return BL_node->parent == nullptr;
    }
    double BL_Node::value() 
    {
        if (points != nullptr){}
        else
        {
            return numeric_limits<double>::infinity();      
        }
        if (!is_leaf()) 
        {
            Pointxy point1 = (*points)[indices.first], point2 = (*points)[indices.second];
            vector<Pointxy> ips = findIntersectionPoints(point1, point2, *sweepline);
            if (ips.size() == 2) 
            {
                if (point1.y < point2.y) 
                {
                    return ips[0].x;
                } 
                else 
                {
                    return ips[1].x;
                }
            } 
            else 
            {
                return ips[0].x;
            }
        } 
        else 
        {
            return (*points)[indices.first].x;
        }
    }


    BL_NPtr find(BL_NPtr root, double x) 
    {
        if (root != nullptr) {}
        else
        {
            return nullptr;
        }
        BL_NPtr BL_node = root;
        while (!BL_node->is_leaf()) 
        {
            if (BL_node->value() >= x) 
            {
                BL_node = BL_node->left;
            } 
            else 
            {
                BL_node = BL_node->right;
            }
        }
        return BL_node;
    }
    int fetch_balance(BL_NPtr BL_node) 
    {
        return find_h(BL_node->left) - find_h(BL_node->right);
    }
    BL_NPtr rotate_left(BL_NPtr BL_node) 
    {
        
        if (BL_node != nullptr){} else
        {
            return nullptr;
        }    
        if (BL_node->right != nullptr){} else
        {
            return BL_node;
        }
        BL_NPtr rnode = BL_node->right;
        if (is_root(BL_node)) {}
        else
        {
            if (BL_node->parent->left != BL_node) 
            {
                BL_node->parent->right = rnode;
                
            } 
            else 
            {
                BL_node->parent->left = rnode;
            }  
        }
        rnode->parent = BL_node->parent;
        BL_node->right = rnode->left;
        if (rnode->left == nullptr) {} else
        {
            rnode->left->parent = BL_node;        
        }
        rnode->left = BL_node;
        BL_node->parent = rnode;
        change_h(BL_node);
        change_h(rnode);
        change_h(rnode->parent);
        return rnode;
    }
    // finds the height
    int find_h(BL_NPtr BL_node) 
    {
        if (BL_node != nullptr)
        {
            return BL_node->height;
        }
        else
        {
            return 0;
        }
    }
    // updates the height
    void change_h(BL_NPtr BL_node) 
    {
        if (BL_node != nullptr)
        {
            BL_node->height = max(find_h(BL_node->left), find_h(BL_node->right)) + 1;
        }
        else
        {
           return;
        }    
    }
    // rotate the tree right around BL_node
    BL_NPtr rotate_right(BL_NPtr BL_node) 
    {
        
        if (BL_node != nullptr){} else
        {
            return nullptr;
        }
        if (BL_node->left != nullptr){} else
        {
            return BL_node;
        }
        BL_NPtr lnode = BL_node->left;
        if (is_root(BL_node)) {}
        else
        {
            if (BL_node->parent->left != BL_node) 
            {
                BL_node->parent->right = lnode;
                
            } 
            else 
            {
                BL_node->parent->left = lnode;
            }   
        }
        lnode->parent = BL_node->parent;
        BL_node->left = lnode->right;
        if (lnode->right == nullptr) {} else
        {
            lnode->right->parent = BL_node;
        }
        lnode->right = BL_node;
        BL_node->parent = lnode;
        change_h(BL_node);
        change_h(lnode);
        change_h(lnode->parent);
        return lnode;
    }
    BL_NPtr make_subtree(int index, int index_behind, double *sweepline,const vector<Pointxy> *points, vector<H_EdgePtr> &edges) 
    {
        BL_NPtr BL_node1 = make_shared<BL_Node>(make_pair(index_behind, index), sweepline, points);
        BL_NPtr BL_node2 = make_shared<BL_Node>(make_pair(index, index_behind), sweepline, points);
        BL_NPtr BL_leaf1 = make_shared<BL_Node>(make_pair(index_behind, index_behind), sweepline, points);
        BL_NPtr BL_leaf2 = make_shared<BL_Node>(make_pair(index, index), sweepline, points);
        BL_NPtr BL_leaf3 = make_shared<BL_Node>(make_pair(index_behind, index_behind), sweepline, points);
        BL_node1->right = BL_node2;
        BL_node2->parent = BL_node1;
        BL_node1->left = BL_leaf1;
        BL_leaf1->parent = BL_node1;
        BL_node2->left = BL_leaf2;
        BL_leaf2->parent = BL_node2;
        BL_node2->right = BL_leaf3;
        BL_leaf3->parent = BL_node2;
        pair<H_EdgePtr, H_EdgePtr> twin_edges = make_twins(index_behind, index);
        BL_node1->edge = twin_edges.first;
        BL_node2->edge = twin_edges.second;
        edges.push_back(twin_edges.first);
        edges.push_back(twin_edges.second);
        connect(BL_leaf1, BL_leaf2);
        connect(BL_leaf2, BL_leaf3);
        change_h(BL_node2);
        change_h(BL_node1);
        return BL_node1;
    }
    // removes a leaf node
    BL_NPtr remove(BL_NPtr leaf) 
    {   
        if (leaf != nullptr){} else
        {
            return nullptr;
        }
        BL_NPtr parent = leaf->parent, grandparent = parent->parent;
        pair<int,int> bp1(leaf->prev->get_id(), leaf->get_id());
        pair<int,int> bp2(leaf->get_id(), leaf->next->get_id());
        pair<int,int> other_bp;
        assert(leaf->next != nullptr);
        assert(leaf->prev != nullptr);
        assert(parent != nullptr);
        assert(grandparent != nullptr);
        assert(parent->has_indices(bp1) || parent->has_indices(bp2));
        if (parent->has_indices(bp1)) 
        {
            other_bp = bp2;
        } 
        else if (parent->has_indices(bp2)) 
        {
            other_bp = bp1;
        }
        BL_NPtr other_subtree;
        if (parent->left != leaf)
        {
            other_subtree = parent->left;
        }
        else
        {
            other_subtree = parent->right;
        }
        other_subtree->parent = grandparent;
        if (grandparent->left != parent) 
        {
            grandparent->right = other_subtree;
        } 
        else 
        {
            grandparent->left = other_subtree;
        }
        BL_NPtr new_root = grandparent;
        // balance the tree
        while (grandparent != nullptr) 
        {
            if (grandparent->has_indices(other_bp))
            {
                grandparent->indices = make_pair(leaf->prev->get_id(), leaf->next->get_id());
            }
            change_h(grandparent);
            int balance = fetch_balance(grandparent);
            if (balance > 1) 
            {
                if (grandparent->left != nullptr && !grandparent->left->is_leaf() && fetch_balance(grandparent->left) < 0) 
                {
                    grandparent->left = rotate_left(grandparent->left);
                }
                grandparent = rotate_right(grandparent);
            } 
            else if (balance < -1) 
            {
                if (grandparent->right != nullptr && !grandparent->right->is_leaf() && fetch_balance(grandparent->right) > 0) 
                {
                    grandparent->right = rotate_right(grandparent->right);
                }
                grandparent = rotate_left(grandparent);
            }
            new_root = grandparent;
            grandparent = grandparent->parent;
        }
        connect(leaf->prev, leaf->next);
        return new_root;
    }

    BL_NPtr replace(BL_NPtr BL_node, BL_NPtr new_node) 
    {   
        if (BL_node != nullptr) {} else
        {
            return new_node;
        }
        double x = new_node->value();
        BL_NPtr parent_node = BL_node->parent;
        new_node->parent = parent_node;
        if (parent_node != nullptr) 
        {
            if (parent_node->value() >= x) 
            {
                parent_node->left = new_node;
            } 
            else 
            {
                parent_node->right = new_node;
            }
        }
        BL_node = new_node;
        while (parent_node != nullptr) 
        {
            change_h(parent_node);
            int balance = fetch_balance(parent_node);
            if (!(balance <= 1)) 
            {
                if (parent_node->left != nullptr && !parent_node->left->is_leaf() && fetch_balance(parent_node->left) < 0) 
                {
                    parent_node->left = rotate_left(parent_node->left);
                }
                parent_node = rotate_right(parent_node);
            } 
            else if (!(balance >= -1)) 
            { 
                if (parent_node->right != nullptr && !parent_node->right->is_leaf() && fetch_balance(parent_node->right) > 0) 
                {
                    parent_node->right = rotate_right(parent_node->right);
                }
                parent_node = rotate_left(parent_node);
            }
            BL_node = parent_node;
            parent_node = parent_node->parent;
        }
        return BL_node;
    }
    BL_NPtr make_simple_subtree(int index, int index_behind, double *sweepline, const vector<Pointxy> *points, vector<H_EdgePtr> &edges) 
    {   
        BL_NPtr BL_node, BL_leaf_l, BL_leaf_r;
        pair<H_EdgePtr, H_EdgePtr> twin_edges = make_twins(index_behind, index);
        edges.push_back(twin_edges.first);
        edges.push_back(twin_edges.second);
        if ((*points)[index].x >= (*points)[index_behind].x) 
        {
            BL_node = make_shared<BL_Node>(make_pair(index_behind, index), sweepline, points);
            BL_leaf_l = make_shared<BL_Node>(make_pair(index_behind, index_behind), sweepline, points);
            BL_leaf_r = make_shared<BL_Node>(make_pair(index, index), sweepline, points);
            BL_node->edge = twin_edges.first;
        } 
        else 
        {
            BL_node = make_shared<BL_Node>(make_pair(index, index_behind), sweepline, points);
            BL_leaf_l = make_shared<BL_Node>(make_pair(index, index), sweepline, points);
            BL_leaf_r = make_shared<BL_Node>(make_pair(index_behind, index_behind), sweepline, points);
            BL_node->edge = twin_edges.second;//twin_edges.first;
        }
        BL_node->left = BL_leaf_l;
        BL_node->right = BL_leaf_r;
        BL_leaf_l->parent = BL_node;
        BL_leaf_r->parent = BL_node;
        connect(BL_leaf_l, BL_leaf_r);
        change_h(BL_node);
        return BL_node;
    }
    pair<BL_NPtr, BL_NPtr> breakpoints(BL_NPtr leaf)
    {   
        if (leaf != nullptr && leaf->next != nullptr && leaf->prev != nullptr){}
        else
        {
            return make_pair<BL_NPtr>(nullptr, nullptr);
        }   
        BL_NPtr parent = leaf->parent, gparent = leaf->parent;
        pair<int,int> bp1(leaf->prev->get_id(), leaf->get_id());
        pair<int,int> bp2(leaf->get_id(), leaf->next->get_id());
        pair<int,int> other_bp;
        bool left_is_missing = true;
        if (parent->has_indices(bp1)) 
        {
            other_bp = bp2;
            left_is_missing = false;
        } 
        else if (parent->has_indices(bp2)) 
        {
            other_bp = bp1;
            left_is_missing = true;
        }
        while (gparent != nullptr) 
        {
            if (!(gparent->has_indices(other_bp))) {}
            else
            {
                break;
            }
            gparent = gparent->parent;
        }
        if (!left_is_missing) 
        {
            return make_pair(parent, gparent);
        } 
        else 
        {
            return make_pair(gparent, parent);
        }
    }


    bool _validate(BL_NPtr BL_node) 
    {   
        if (BL_node != nullptr){} else
        {
            return true;
        }
        if (!BL_node->is_leaf()) 
        {
            if (BL_node->left == nullptr || BL_node->right == nullptr) 
            {
                cout << " BP WITHOUT LEAF: " << BL_node->indices.first << ", " << BL_node->indices.second << endl;
                return false;
            }
        } 
        else 
        {
            if (BL_node->left != nullptr || BL_node->right != nullptr) 
            {
                cout << "LEAF NOT A LEAF: " << BL_node->indices.first << ", " << BL_node->indices.second << endl;
                return false;
            }
        }
        return true;
    }

    bool _check_balance(BL_NPtr BL_node) 
    {
        if (BL_node != nullptr){}
        else
        {
            return true;
        }
        if (_check_balance(BL_node->left) && _check_balance(BL_node->right)) 
        {
            if (fabs(fetch_balance(BL_node)) <= 1){}
            else
            {
                cout << "+unbalanced (" << BL_node->indices.first << ", " << BL_node->indices.second << ")" << endl;
                return false;      
            }
        }
        return true;
    }
}

// Functions for DCEL
namespace DCEL {  
    Vertex::Vertex(const Pointxy &pos, H_EdgePtr incident_edge) : point(pos), edge(incident_edge)
    {}
    H_edge::H_edge(int _l_index, int _r_index, VertexPtr _vertex) :
        l_index(_l_index), r_index(_r_index), vertex(_vertex) {}
    H_EdgePtr H_edge::vertexNextCCW() 
    {
        return twin->prev;
    }
    H_EdgePtr H_edge::vertexNextCW() 
    {
        return next->twin;
    }
    pair<H_EdgePtr, H_EdgePtr> make_twins(int left_index, int right_index) 
    {   
        H_EdgePtr h_ptr = make_shared<H_edge>(left_index, right_index);
        H_EdgePtr h_ptr_twin = make_shared<H_edge>(right_index, left_index);
        
        h_ptr->twin = h_ptr_twin;
        h_ptr_twin->twin = h_ptr;
        
        return make_pair(h_ptr, h_ptr_twin);
    }
    pair<H_EdgePtr, H_EdgePtr> make_twins(const pair<int,int> &indices) 
    {
        
        return make_twins(indices.first, indices.second);
    }
    void connect_h_edges(H_EdgePtr point1, H_EdgePtr point2) 
    {
        point1->next = point2;
        point2->prev = point1;
    }
}


// Functions for circle
bool findCircleCenter(const Pointxy &point1, const Pointxy &point2, const Pointxy &point3, Pointxy &center) 
{   
    Pointxy p_u1 = (point1 - point2).normalized(), p_u2 = (point3 - point2).normalized();   
    double cross_prod = vector_prod(p_u1, p_u2);
    if (fabs(cross_prod) >= CIRCLE_CENTER_EPSILON) {} else
    {
        return false;
    }    
    Pointxy pc1=0.5*(point1+point2), pc2=0.5*(point2+point3);
    double scal_b1 = scalar_prod(p_u1, pc1), scal_b2 = scalar_prod(p_u2, pc2);
    center.x = (scal_b1*p_u2.y-scal_b2*p_u1.y)/cross_prod;
    center.y = (p_u1.x*scal_b2-p_u2.x*scal_b1)/cross_prod;
    return true;
}

struct Event 
{
    enum { SITE = 0, CIRCLE = 1, SKIP = 2, };
    int type;
    Pointxy point;
    int index;
    Pointxy center;
    bl::BL_NPtr arc;
    Event(int _index = -1, int _type = Event::SKIP, const Pointxy &_point = Pointxy(0.0, 0.0)) :
    index(_index), type(_type), point(_point), arc(nullptr) {}
};
// function definitions for voronai
struct Point2DComparator 
{
    bool operator()(const Pointxy &point1, const Pointxy &point2) 
    {
        return (point1.y == point2.y && point1.x > point2.x) || point1.y > point2.y;
    }
};
shared_ptr<Event> checkCircleEvent(bl::BL_NPtr n1, bl::BL_NPtr n2, bl::BL_NPtr n3,const vector<Pointxy> &points, double sweepline) 
{  
    if (n1 != nullptr && n2 != nullptr && n3 != nullptr){}
    else
    {
        return nullptr;
    }   
    Pointxy point1 = points[n1->get_id()];
    Pointxy point2 = points[n2->get_id()];
    Pointxy point3 = points[n3->get_id()];
    Pointxy center, bottom;
    if (point2.y <= point1.y || point2.y <= point3.y){}
    else
    {
        return nullptr;
    }
    if (findCircleCenter(point1, point2, point3, center)){}
    else
    {
        return nullptr;
    }
    bottom = center;
    bottom.y += (center - point2).norm();
    if (fabs(bottom.y - sweepline) < POINT_EPSILON || sweepline < bottom.y) 
    {
        shared_ptr<Event> e = make_shared<Event>(-1, Event::CIRCLE, bottom);
        e->center = center;
        e->arc = n2;
        n2->circle_event = e;
        return e;
    }
    return nullptr;
}
struct Point2DComparator2 {
    bool operator()(const Pointxy &point1, const Pointxy &point2) 
    {
        return (point1.y == point2.y && point1.x < point2.x) || point1.y < point2.y;
    }
};
struct EventPtrComparator 
{
    Point2DComparator point_cmp;
    bool operator()(const shared_ptr<Event> &e1, const shared_ptr<Event> &e2) 
    {
        return point_cmp(e1->point, e2->point);
    }
};

void build_voronoi(const vector<Pointxy> &points, vector<bl::H_EdgePtr> &halfedges, vector<bl::VertexPtr> &vertices, vector<bl::H_EdgePtr> &faces) {
    priority_queue<shared_ptr<Event>, vector<shared_ptr<Event>>, EventPtrComparator> pq;
    for (size_t i = 0; i < points.size(); ++i) 
    {
        pq.push(make_shared<Event>(static_cast<int>(i), Event::SITE, points[i]));
    }
    faces.resize(points.size(), nullptr);
    bl::BL_NPtr root;
    double sweepline = 0L;
    while (!pq.empty()) 
    {
        shared_ptr<Event> e = pq.top(); pq.pop();
        sweepline = e->point.y;
        if (e->type == Event::SITE) 
        {
            int point_i = e->index;
            if (root != nullptr) 
            {
                bl::BL_NPtr arc = bl::find(root, e->point.x);
                bl::BL_NPtr subtree, left_leaf, right_leaf;   
                if (arc->circle_event != nullptr) 
                {
                    shared_ptr<Event> circle_e = arc->circle_event;
                    circle_e->type = Event::SKIP;
                }
                int isp_num = intersectionPoints(points[arc->get_id()], e->point, sweepline);
                if (isp_num == 1) 
                {
                    subtree = bl::make_simple_subtree(point_i, arc->get_id(), &sweepline, &points, halfedges);
                    left_leaf = subtree->left;
                    right_leaf = subtree->right;
                } 
                else if (isp_num == 2) 
                {
                    subtree = bl::make_subtree(point_i, arc->get_id(), &sweepline, &points, halfedges);
                    left_leaf = subtree->left;
                    right_leaf = subtree->right->right;
                } 
                else 
                {
                    continue;
                }
                if (arc->prev != nullptr)
                {
                    bl::connect(arc->prev, left_leaf);
                }
                if (arc->next != nullptr)
                {
                    bl::connect(right_leaf, arc->next);
                }
                root = bl::replace(arc, subtree);
                shared_ptr<Event> circle_event = checkCircleEvent(left_leaf->prev, left_leaf, left_leaf->next, points, sweepline);
                if (circle_event != nullptr) 
                {
                    pq.push(circle_event);
                }
                circle_event = checkCircleEvent(right_leaf->prev, right_leaf, right_leaf->next, points, sweepline);
                if (circle_event != nullptr) 
                {
                    pq.push(circle_event);
                }
            } 
            else 
            {
                root = make_shared<bl::BL_Node>(make_pair(point_i, point_i), &sweepline, &points);
            }
        } 
        else if (e->type == Event::CIRCLE) 
        {   
            bl::BL_NPtr arc = e->arc, prev_leaf, next_leaf;
            pair<bl::BL_NPtr, bl::BL_NPtr> breakpoints = bl::breakpoints(arc);
            if (breakpoints.first != nullptr && breakpoints.second != nullptr) {}
            else
            {
                continue;
            }
            double v1 = breakpoints.first->value(), v2 = breakpoints.second->value();
            if (fabs(v1 - v2) <= BREAKPOINTS_EPSILON){}
            else
            {
                continue;
            }
            bl::VertexPtr vertex = make_shared<bl::Vertex>(e->center);
            bl::H_EdgePtr h_first = breakpoints.first->edge;
            bl::H_EdgePtr h_second = breakpoints.second->edge;
            vertices.push_back(vertex);
            if (arc->prev != nullptr && arc->prev->circle_event != nullptr) 
            {
                shared_ptr<Event> circle_e = arc->prev->circle_event;
                circle_e->type = Event::SKIP; // ignore corresponding event
            }
            if (arc->next != nullptr && arc->next->circle_event != nullptr) 
            {
                shared_ptr<Event> circle_e = arc->next->circle_event;
                circle_e->type = Event::SKIP; // ignore corresponding event
            }
            prev_leaf = arc->prev;
            next_leaf = arc->next;
            assert(prev_leaf != nullptr);
            assert(next_leaf != nullptr);
            bl::BL_NPtr new_edge_node;
            if (arc->parent != breakpoints.first)
            {
                new_edge_node = breakpoints.first; 
            }
            else
            {
                new_edge_node = breakpoints.second;
            }
            root = bl::remove(arc);
            pair<bl::H_EdgePtr, bl::H_EdgePtr> twin_nodes = bl::make_twins(prev_leaf->get_id(), next_leaf->get_id());
            new_edge_node->edge = twin_nodes.first;
            bl::connect_h_edges(h_second, h_first->twin);
            bl::connect_h_edges(h_first, twin_nodes.first);
            bl::connect_h_edges(twin_nodes.second, h_second->twin);
            h_first->vertex = vertex;
            h_second->vertex = vertex;
            twin_nodes.second->vertex = vertex;
            vertex->edge = h_second;
            halfedges.push_back(twin_nodes.first);
            halfedges.push_back(twin_nodes.second);
            if (prev_leaf != nullptr && next_leaf != nullptr) 
            {
                shared_ptr<Event> circle_event = checkCircleEvent(prev_leaf->prev, prev_leaf, next_leaf, points, sweepline);
                if (circle_event == nullptr){}
                else
                {
                    pq.push(circle_event);
                }
                circle_event = checkCircleEvent(prev_leaf, next_leaf, next_leaf->next, points, sweepline);
                if (circle_event == nullptr){}
                else
                {
                    pq.push(circle_event);
                }
            }
        }
    }
    for (size_t i = 0; i < halfedges.size(); ++i) 
    {
        bl::H_EdgePtr he = halfedges[i];
        if (he->prev == nullptr || faces[he->l_index] == nullptr) 
        {
            faces[he->l_index] = he;
        }
    }
}
