
// Name:Akash Tiwari
// Roll no:17CS10003
// NOTE: AT N>25, the lines in .svg become very light.

#include "Fortune_header.hpp"
#include<bits/stdc++.h>
#include<memory>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#define MAX 300

double x_max=0,y_max=0;
double x_min=0,y_min=0;
namespace bl = beachline;

using namespace std;


// draws a circle in svg
void draw_circle(const Pointxy &c, double r, ofstream & fout , double shift_x, double shift_y) 
{
    fout<<"<circle cx=\""<<c.x+shift_x<<"\" "<<"cy=\""<<c.y+shift_y<<"\" "<<"r=\""<<r<<"\" stroke=\"#66ff33\" stroke-width=\"0\" fill=\"#66ff33\" fill-opacity=\"0.4\" />\n";   
}

void pin_point(const Pointxy &c, ofstream & fout , double shift_x, double shift_y)
{       
 fout<<"<circle cx=\""<<c.x+shift_x<<"\" "<<"cy=\""<<c.y+shift_y<<"\" "<<"r=\""<<2<<"\" stroke=\"black\" stroke-width=\"0\" fill=\"#000000\" />\n";
}

// random point(x and y co-ords are doubles) generator
vector<Pointxy> randomPoint(int num) 
{
	// providing seed as time for rand
    srand(static_cast<unsigned int>(time(0)));
    // list of points to return
    vector<Pointxy> P_list;
    for (int itr = 0; itr < num; ++itr) 
    {
    	//Finding random co-ordinates
        double x = (double(rand())*MAX)/ double(RAND_MAX);
        double y = (double(rand())*MAX)/ double(RAND_MAX);
        int id=itr+1;
        // defining point to be pushed
        Pointxy p_rand;
        p_rand.x=x;
        p_rand.y=y;
        p_rand.id=id;
        P_list.push_back(p_rand);
    }
    return P_list;
}

void InitEdgeVPoints(bl::H_EdgePtr h, vector<double> &x,vector<double> &y,const vector<Pointxy> &points) 
{ 
    if (h->vertex != nullptr && h->twin->vertex != nullptr) 
    {
        x[0] = h->vertex->point.x;
        x[1] = h->twin->vertex->point.x;
        y[0] = h->vertex->point.y;
        y[1] = h->twin->vertex->point.y;
      auto  temp_max=max_element(x.begin(),x.end());
      auto  temp_min=min_element(x.begin(),x.end());        
        if(x_max< *temp_max)
        {
            x_max= (*temp_max);
        }
        if(x_min> *temp_min)
        {
            x_min= (*temp_min);
        }
        temp_max=max_element(y.begin(),y.end());
        temp_min=min_element(y.begin(),y.end());        
        if(y_max< *temp_max)
        {
            y_max= (*temp_max);
        }
        if(y_min> *temp_min)
        {
            y_min= (*temp_min);         
        }   
    } else if (h->vertex != nullptr) 
    { 
        x[0] = h->vertex->point.x;
        y[0] = h->vertex->point.y;
        Pointxy normal = (points[h->l_index] - points[h->r_index]).normalized().Rotate90CCW();
        x[1] = x[0] + normal.x * 1000;
        y[1] = y[0] + normal.y * 1000;
      auto  temp_max=max_element(x.begin(),x.end());
      auto  temp_min=min_element(x.begin(),x.end());        
        if(x_max< *temp_max)
        {
            x_max= (*temp_max);
        }
        if(x_min> *temp_min)
        {
            x_min= (*temp_min);
        }
        temp_max=max_element(y.begin(),y.end());
        temp_min=min_element(y.begin(),y.end());        
        if(y_max< *temp_max)
         {
            y_max= (*temp_max);
         }
        if(y_min> *temp_min)
        {
            y_min= (*temp_min);         
        }  
    } else if (h->twin->vertex != nullptr) 
    {     
        x[0] = h->twin->vertex->point.x;
        y[0] = h->twin->vertex->point.y;
        Pointxy normal = (points[h->twin->l_index] - points[h->twin->r_index]).normalized().Rotate90CCW();
        x[1] = x[0] + normal.x * 1000;
        y[1] = y[0] + normal.y * 1000;
      auto  temp_max=max_element(x.begin(),x.end());
      auto  temp_min=min_element(x.begin(),x.end());        
        if(x_max< *temp_max)
         {
            x_max= (*temp_max);
         }
        if(x_min> *temp_min)
         {
            x_min= (*temp_min);
         }
        temp_max=max_element(y.begin(),y.end());
        temp_min=min_element(y.begin(),y.end());        
        if(y_max< *temp_max)
         {
            y_max= (*temp_max);
         }
        if(y_min> *temp_min)
         {
            y_min= (*temp_min); 
         }  
    } 
    else 
    {   
        Pointxy p1 = points[h->l_index], p2 = points[h->r_index];  
        Pointxy normal = (p1 - p2).normalized().Rotate90CCW();
        Pointxy c = 0.5 * (p1 + p2); 
        x[0] = c.x + normal.x *1000 ;
        x[1] = c.x - normal.x *1000 ;
        y[0] = c.y + normal.y * 1000;
        y[1] = c.y - normal.y * 1000;
     auto   temp_max=max_element(x.begin(),x.end());
     auto   temp_min=min_element(x.begin(),x.end());        
        if(x_max< *temp_max)
         {
            x_max= (*temp_max);
         }
        if(x_min> *temp_min)
         {
            x_min= (*temp_min);
         }
        temp_max=max_element(y.begin(),y.end());
        temp_min=min_element(y.begin(),y.end());        
        if(y_max< *temp_max)
        {
           y_max= (*temp_max);
        }
        if(y_min> *temp_min)
        {
            y_min= (*temp_min); 
        }
    }
}


int main() {
    int n;
    cout<<"Enter the number of sites in the diagram:\n";
    cin>>n;
    // generate a list of random points
    vector<Pointxy> points = randomPoint(n);
   vector<bl::H_EdgePtr> h_edges, faces;
   vector<bl::VertexPtr> vertices;
    ofstream fout;
    fout.open("t4.svg");
   build_voronoi(points, h_edges, vertices, faces);
    if(n<=1)
     {
       cout<<"Cannot draw for n<=1\n";
       return 0;
     }
     else if(n==2)
     {
        for (size_t itr = 0; itr < h_edges.size(); ++itr) 
        {
          bl::H_EdgePtr h = h_edges[itr];
          vector<double> x(2, 0.0), y(2, 0.0);

        InitEdgeVPoints(h, x, y, points);
        }
  double x_min_1=x_min;
  double x_max_1=x_max;
  double y_min_1=y_min;
  double y_max_1=y_max;  
  double shift_x=x_max_1-x_min_1;
  double shift_y=y_max_1-y_min_1;
   fout<<"<svg xmlns=\"http://www.w3.org/2000/svg\">\n<rect width=\""<<2*shift_x+100<<"\" height=\""<<2*shift_y+100<<"\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />\n";
   for(size_t i = 0; i < points.size(); ++i)
    {
     pin_point(points[i],fout,shift_x,shift_y);
    }
      for (size_t i = 0; i < h_edges.size(); ++i) 
      {
        bl::H_EdgePtr h = h_edges[i];
        vector<double> x(2, 0.0), y(2, 0.0);
        InitEdgeVPoints(h, x, y, points);
        fout<<"<line x1=\""<<x[0]+shift_x<<"\" y1=\""<<y[0]+shift_y<<"\" x2=\""<<x[1]+shift_x<<"\" y2=\""<<y[1]+shift_y<<"\" style=\"stroke:rgb(255,0,0);stroke-width:1\" />\n";
      } 
      double dist2;
      Pointxy diff = points[0]-points[1];
      dist2= diff.norm();
      // draw the 2 circles
      draw_circle(points[0],dist2/2.0,fout,shift_x,shift_y);
      draw_circle(points[1],dist2/2.0,fout,shift_x,shift_y);
      fout<<"</svg>\n"; 
     }
    else
    {
   for (size_t i = 0; i < h_edges.size(); ++i) 
   {
        bl::H_EdgePtr h = h_edges[i];
        vector<double> x(2, 0.0), y(2, 0.0);
        InitEdgeVPoints(h, x, y, points);
    }
  double x_min_1=x_min;
  double x_max_1=x_max;
  double y_min_1=y_min;
  double y_max_1=y_max;  
  double shift_x=x_max_1-x_min_1;
  double shift_y=y_max_1-y_min_1; 
   fout<<"<svg xmlns=\"http://www.w3.org/2000/svg\">\n<rect width=\""<<2*shift_x+100<<"\" height=\""<<2*shift_y+100<<"\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />\n";
   for(size_t i = 0; i < points.size(); ++i)
    {
     pin_point(points[i],fout,shift_x,shift_y);
    }
      for (size_t i = 0; i < h_edges.size(); ++i) 
      {
       bl::H_EdgePtr h = h_edges[i];
       vector<double> x(2, 0.0), y(2, 0.0);
       InitEdgeVPoints(h, x, y, points);
       fout<<"<line x1=\""<<x[0]+shift_x<<"\" y1=\""<<y[0]+shift_y<<"\" x2=\""<<x[1]+shift_x<<"\" y2=\""<<y[1]+shift_y<<"\" style=\"stroke:rgb(255,0,0);stroke-width:1\" />\n";
      }

     // draws the circle with least radius for each site 
    for(size_t itr = 0; itr < faces.size(); ++itr) 
    { 
      Pointxy site =points[itr];
      bl::H_EdgePtr h = faces[itr];
      bl::H_EdgePtr h1=h;
      while(h1->vertex!=nullptr)
      {h1=h1->next;
        if(h1==h)
        {
          break;
        }
      } 
      double min_d=-1;
      double dist;
        do 
        {
          Pointxy site_2 = points[h1->r_index];
          Pointxy diff = site-site_2;
            dist= diff.norm();
            if(min_d>dist||min_d==-1)
            {
              min_d = dist;
            }
            h1=h1->prev;
        }while(h1->twin->vertex!=nullptr&&h1!=h);
        Pointxy site_2 = points[h1->r_index];
         Pointxy diff = site-site_2;
            dist= diff.norm();
            if(min_d>dist||min_d==-1)
            {
              min_d = dist;
            }        
      draw_circle(site,min_d/2.0,fout,shift_x,shift_y);
    }
    
    cout<<"t4.svg has been generated.\n";
    fout<<"</svg>\n"; 
    return 0;
    }
}

