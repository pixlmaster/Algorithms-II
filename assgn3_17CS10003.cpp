#include <iostream>
#include <bits/stdc++.h>

using namespace std;

struct point{
	int x;
	int y;
};

int comparex(point a, point b){
	if (a.x<=b.x){
		return true;
	}
	else{
		return false;
	}
}

int comparey(point a, point b){
	if (a.y<=b.y){
		return true;
	}
	else{
		return false;
	}
}

float dist(point p1, point p2) 
{ 
	return sqrt( (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)); 
} 

float brute(point p[], int n) 
{ 
	float min = 999999999; 
	for (int i = 0; i < n; ++i) 
		for (int j = i+1; j < n; ++j) 
			if (dist(p[i], p[j]) < min) 
				min = dist(p[i], p[j]); 
	return min; 
} 

float closstrip(point strip[], int s, float d) 
{ 
	float min = d; 

	for (int i = 0; i < s; ++i) 
		for (int j = i+1; j < s && (strip[j].y - strip[i].y) < min; ++j) 
			if (dist(strip[i],strip[j]) < min) 
				min = dist(strip[i], strip[j]); 

	return min; 
} 

float closerec(point sortx[], point sorty[] , int n){
	if(n<4){
		return brute(sortx,n);
	}
	int middle = n/2;
	point pmiddle = sortx[middle];

	point sortyl[middle+1];
	point sortyr[n-middle-1];
	int indl=0, indr=0;
	for (int i = 0; i < n; ++i)
	{
		if (sorty[i].x<=pmiddle.x)
		{
			sortyl[indl] = sorty[i];
			indl++;
		}
		else{
			sortyr[indr] = sorty[i];
			indr++;
		}
	}

	float dl = closerec(sortx,sortyl,middle);
	float dr = closerec(sortx+middle , sortyr , n-middle);
	float d;
	if(dl<dr){
		d=dl;
	}
	else{
		d=dr;
	}
	point strip[n];
	int j = 0;
	for (int i = 0; i < n; ++i)
	{
		if (abs(sorty[i].x - pmiddle.x) < d)
		{
			strip[j] = sorty[i];
			j++;
		}
	}
	float stripd = closstrip(strip,j,d);
	if(dl<stripd){
		return d;
	}
	else{
		return stripd;
	}
}

int main(){
	int n;
	cin>>n;
	point p[n];
	point sortx[n];
	point sorty[n];
	srand(time(0));
	for (int i = 0; i < n; ++i)
	{
		int randx = rand()%1000;
		int randy = rand()%750;
		p[i].x=randx;
		p[i].y=randy;
		sortx[i]=p[i];
		sorty[i]=p[i];
	}
	sort(sortx,sortx+n, comparex);
	sort(sorty,sorty+n, comparey);
	float radius = closerec(sortx,sorty,n);
	ofstream fout;
	fout.open("ans.svg");
	fout<<"<svg xmlns=\"http://www.w3.org/2000/svg\">"<<endl;
	fout<<"<rect width=\"2000\" height=\"1500\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />"<<endl;
	for (int i = 0; i < n; ++i)
	{
		fout<<"<circle cx=\""<<p[i].x<<"\" cy=\""<<p[i].y<<"\" r=\""<<radius/2<<"\" stroke=\"black\" stroke-width=\"0\" fill=\"#000000\" fill-opacity=\"0.8\" />"<<endl;
		fout<<"<text x=\""<<p[i].x<<"\" y=\""<<p[i].y<<"\">"<<i+1<<"</text>"<<endl;
	}
	fout<<"</svg>"<<endl;
	fout.close();
}