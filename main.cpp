#if defined(UNICODE) && !defined(_UNICODE)
#define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
#define UNICODE
#endif

#include <tchar.h>
#include <windows.h>
#include<bits/stdc++.h>
#include "math.h"
#define MAXENTRIES 600

using namespace std;
LPCSTR mouse=IDC_ARROW;
HBRUSH hBrush = CreateSolidBrush(RGB(255, 255, 255)); // RGB(255, 255, 255) represents white color


/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure (HWND, UINT, WPARAM, LPARAM);
void load(HWND, HDC &);
void save(HWND &);

/*  Make the class name into a global variable  */
TCHAR szClassName[ ] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain (HINSTANCE hThisInstance,
                    HINSTANCE hPrevInstance,
                    LPSTR lpszArgument,
                    int nCmdShow)
{
    HWND hwnd;               /* This is the handle for our window */
    MSG messages;            /* Here messages to the application are saved */
    WNDCLASSEX wincl;        /* Data structure for the windowclass */

    /* The Window structure */
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = szClassName;
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof (WNDCLASSEX);

    /* Use default icon and mouse-pointer */
    wincl.hIcon = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hIconSm = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hCursor = LoadCursor (NULL, IDC_ARROW);
    wincl.lpszMenuName = NULL;                 /* No menu */
    wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
    wincl.cbWndExtra = 0;                      /* structure or the window instance */
    /* Use Windows's default colour as the background of the window */
    wincl.hbrBackground = hBrush;

    /* Register the window class, and if it fails quit the program */
    if (!RegisterClassEx (&wincl))
        return 0;

    /* The class is registered, let's create the program*/
    hwnd = CreateWindowEx (
               0,                   /* Extended possibilites for variation */
               szClassName,         /* Classname */
               _T("Code::Blocks Template Windows App"),       /* Title Text */
               WS_OVERLAPPEDWINDOW, /* default window */
               CW_USEDEFAULT,       /* Windows decides the position */
               CW_USEDEFAULT,       /* where the window ends up on the screen */
               800,                 /* The programs width */
               600,                 /* and height in pixels */
               HWND_DESKTOP,        /* The window is a child-window to desktop */
               NULL,                /* No menu */
               hThisInstance,       /* Program Instance handler */
               NULL                 /* No Window Creation data */
           );

    /* Make the window visible on the screen */
    ShowWindow (hwnd, nCmdShow);

    /* Run the message loop. It will run until GetMessage() returns 0 */
    while (GetMessage (&messages, NULL, 0, 0))
    {
        /* Translate virtual-key messages into character messages */
        TranslateMessage(&messages);
        /* Send message to WindowProcedure */
        DispatchMessage(&messages);
    }

    /* The program return-value is 0 - The value that PostQuitMessage() gave */
    return messages.wParam;
}

/*______________________________________ Line _________________________________________*/

void swap(int &x, int &y)
{
    int temp = x;
    x = y;
    y = temp;
}

int max(int x, int y)
{
    return (x>y? x: y);
}

void lineDDA(HDC hdc,int x1,int y1,int x2, int y2, COLORREF c)
{
    int dx = x2-x1, dy = y2-y1;
    double m = (double)dy / dx,mInv = (double)dx / dy;
    if(abs(dx) >= abs(dy))
    {
        if(x1>x2)
        {
            swap(x1,x2);
            swap(y1,y2);
        }
        double y = y1;
        int x = x1;
        SetPixel(hdc,x1,y1,c);
        while(x < x2)
        {
            y += m;
            x++;
            SetPixel(hdc,x,round(y),c);
        }
    }
    else
    {
        if(y1>y2)
        {
            swap(x1,x2);
            swap(y1,y2);
        }
        int y = y1;
        double x = x1;
        SetPixel(hdc,x1,y1,c);
        while(y < y2)
        {
            y ++;
            x += mInv;
            SetPixel(hdc,round(x),y,c);
        }
    }
}

void lineBresenham(HDC hdc,int x1,int y1,int x2, int y2, COLORREF c)
{
    int dx = x2 - x1, dy = y2 - y1;
    int x, y;
    int d, change1, change2;
    if(abs(dx) >= abs(dy))
    {
        if(x1>x2)
        {
            swap(x1,x2);
            swap(y1,y2);
        }
        dx = x2 - x1;
        dy = y2 - y1;
        x = x1 ;
        y = y1 ;
        d = dx - 2 * abs(dy);   // d initial value at (x1,y1)
        change1 = -2 * abs(dy);
        change2 = 2 * (dx - abs(dy));
        while(x <= x2)
        {
            SetPixel(hdc,x,y,c);
            if(d >= 0)
            {
                d+=change1;
            }
            else
            {
                // x increases , but y ( maybe yes , maybe no )
                if(dy > 0)
                    y++;
                else
                    y--;
                d+=change2;
            }
            x++;
        }
    }
    else
    {
        if(y1>y2)
        {
            swap(x1,x2);
            swap(y1,y2);
        }
        dx = x2 - x1;
        dy = y2 - y1;
        x = x1 ;
        y = y1 ;
        d = 2 * abs(dx) - dy;   // d initial value at (x1,y1)
        change1 = 2 * abs(dx);
        change2 = 2 * (abs(dx) - dy);
        while(y <= y2)
        {
            SetPixel(hdc,x,y,c);
            if(d > 0)
            {
                d+=change2;
                if(dx > 0)
                    x++;
                else
                    x--;
            }
            else
            {
                d+=change1;
            }
            y++;
        }
    }
}



void ParametricLine(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double t=0;
    while(t<=1)
    {
        int x=x1 + round(dx * t);
        int y = y1 + round(dy * t);
        SetPixel(hdc, x, y, color);
        t+=0.001;
    }

}
/*______________________________________ Circle _________________________________________*/


int distance(int xc,int yc,int xr,int yr)
{
    return (int)sqrt(pow(abs(xr-xc),2)+pow(abs(yr-yc),2)) ;
}
void cicrlepolar(HDC hdc,int xc,int yc,int R,COLORREF c)
{
    double dtheta=1.0/R;
    for(double theta=0; theta<6.28; theta+=dtheta)
    {
        double x=xc+R*cos(theta);
        double y=yc+R*sin(theta);
        SetPixel(hdc,round(x),round(y),c);

    }
}

void Draw8Points(HDC hdc,int xc,int yc, int a, int b,COLORREF color)
{
    SetPixel(hdc, xc+a, yc+b, color);
    SetPixel(hdc, xc-a, yc+b, color);
    SetPixel(hdc, xc-a, yc-b, color);
    SetPixel(hdc, xc+a, yc-b, color);
    SetPixel(hdc, xc+b, yc+a, color);
    SetPixel(hdc, xc-b, yc+a, color);
    SetPixel(hdc, xc-b, yc-a, color);
    SetPixel(hdc, xc+b, yc-a, color);
}
void CircleIterativePolar(HDC hdc,int xc,int yc, int R,COLORREF color)
{
    double dtheta=1.0/R;
    double ct=cos(dtheta),st=sin(dtheta);

    double x=R,y=0;
    Draw8Points(hdc,xc,yc,R,0,color);
    while(x>y)
    {
        double x1=x*ct-y*st;
        y=x*st+y*ct;
        x=x1;
        Draw8Points(hdc,xc,yc,round(x),round(y),color);
    }
}
void CircleDirect(HDC hdc,int xc,int yc, int R,COLORREF color)
{
    int x=0;
    double y=R;
    Draw8Points(hdc,xc,yc,0,R,color);
    while(x<y)
    {
        x++;
        y=sqrt((double)(R*R-x*x));
        Draw8Points(hdc,xc,yc,round(x),round(y),color);

    }


}
void Circlemidpoint(HDC hdc,int xc,int yc, int R,COLORREF color)
{
    int x=0,y=R;
    Draw8Points(hdc,xc,yc,x,y,color);
    while(x<y)
    {
        double d=pow(x+1,2)+pow(y-1/2,2)-(R*R);
        if(d<=0)
        {
            x++;
        }
        else
        {
            x++;
            y--;
        }
        Draw8Points(hdc,xc,yc,x,y,color);


    }
}
void CircleModifiedMidpoint(HDC hdc,int xc,int yc, int R,COLORREF color)
{
    int x=0,y=R, d=1-R;
    Draw8Points(hdc,xc,yc,x,y,color);
    while(x<y)
    {
        if(d<0)
        {
            d+=2*x+3;
            x++;
        }
        else
        {

            d+=2*(x-y)+5;
            x++;
            y--;
        }

        Draw8Points(hdc,xc,yc,x,y,color);
    }
}

/*______________________________________ Fill the shape ________________________________________*/

//Circle with lines

void CircleFillWithLines(HDC hdc, int xc, int yc, int a,int b,COLORREF color, int quarter)
{
    // Draw the border of the circle

    // Draw lines to fill the desired quarter
    switch (quarter)
    {

    case 1:
        lineDDA(hdc, xc, yc, xc + a, yc - b, color);
        lineDDA(hdc, xc, yc, xc + b, yc - a, color);

        break;
    case 2:
        lineDDA(hdc, xc, yc, xc + a, yc + b, color);
        lineDDA(hdc, xc, yc, xc + b, yc + a, color);
        break;
    case 3:
        lineDDA(hdc, xc, yc, xc - a, yc + b, color);
        lineDDA(hdc, xc, yc, xc - b, yc + a, color);
        break;
    case 4:
        lineDDA(hdc, xc, yc, xc - b, yc - a, color);
        lineDDA(hdc, xc, yc, xc - a, yc - b, color);
        break;
    default:
        break;
    }
}
//Fill circle quarter with lines

void fillCircleWithLines(HDC hdc, int xc, int yc, int R,  COLORREF c,int quarter )
{
    int x=0,y=R, d=1-R;
    Draw8Points(hdc,xc,yc,x,y,c);
    while(x<y)
    {
        if(d<0)
        {
            d+=2*x+3;
            x++;
        }
        else
        {

            d+=2*(x-y)+5;
            x++;
            y--;
        }
        Draw8Points(hdc, xc, yc, round(x), round(y), c);
        CircleFillWithLines (hdc, xc, yc, round(x), round(y),  c,quarter);
    }
}
//Fill circle quarter with other circles

void draw8Points(HDC hdc, int xc, int yc, int x, int y, COLORREF c, int quarter)
{
    switch (quarter)
    {
    case 0:
        SetPixel(hdc, xc + x, yc + y, c);
        SetPixel(hdc, xc - x, yc + y, c);
        SetPixel(hdc, xc + x, yc - y, c);
        SetPixel(hdc, xc - x, yc - y, c);

        SetPixel(hdc, xc + y, yc + x, c);
        SetPixel(hdc, xc - y, yc + x, c);
        SetPixel(hdc, xc + y, yc - x, c);
        SetPixel(hdc, xc - y, yc - x, c);
        break;
    case 1:
        SetPixel(hdc, xc + x, yc - y, c);
        SetPixel(hdc, xc + y, yc - x, c);
        break;
    case 2:
        SetPixel(hdc, xc - x, yc - y, c);
        SetPixel(hdc, xc - y, yc - x, c);
        break;
    case 3:
        SetPixel(hdc, xc - x, yc + y, c);
        SetPixel(hdc, xc - y, yc + x, c);
        break;
    case 4:
        SetPixel(hdc, xc + x, yc + y, c);
        SetPixel(hdc, xc + y, yc + x, c);
        break;
    default:
        break;
    }
}

//Fill Circle with Circles

void DrawCircle(HDC hdc, int xc, int yc, int radius, COLORREF color, int quarter)
{
    double dTheta = 1.0/radius, x, y;
    for (double theta = 0; theta < M_PI/2; theta += dTheta)
    {
        x = radius * cos(theta);
        y = radius * sin(theta);
        draw8Points(hdc, xc, yc, x, y, color, quarter);
    }
}

void CircleFillWithCircles(HDC hdc, int xc, int yc, int radius, COLORREF color, int quarter)
{
    for (int r = 0; r <= radius; r++)
    {
        DrawCircle(hdc, xc, yc, r, color, quarter);
    }
    DrawCircle(hdc, xc, yc, radius, color, 0);
}



//Hermite curve
struct Vector
{
    double v[2];
    Vector(double x = 0, double y = 0)
    {
        v[0] = x;
        v[1] = y;
    }
    const double& operator[](int i) const
    {
        return v[i];
    }
};


class Vector4
{
    double v[4];
public:
    Vector4(double a=0,double b=0,double c=0,double d=0)
    {
        v[0]=a;
        v[1]=b;
        v[2]=c;
        v[3]=d;
    }


    Vector4(double a[])
    {
        memcpy(v,a,4*sizeof(double));
    }
    double& operator[](int i)
    {
        return v[i];
    }
};
class Matrix4
{
    Vector4 M[4];
public:
    Matrix4(double A[])
    {
        memcpy(M,A,16*sizeof(double));
    }
    Vector4& operator[](int i)
    {
        return M[i];
    }
};


Vector4 operator*(Matrix4 M,Vector4& b) // right multiplication of M by b
{
    Vector4 res;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            res[i]+=M[i][j]*b[j];

    return res;
}


double DotProduct(Vector4& a,Vector4& b) //multiplying a raw vector by a column vector
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
}

Vector4 GetHermiteCoeff(double x0,double s0,double x1,double s1)
{
    static double H[16]= {2,1,-2,1,-3,-2,3,-1,0,1,0,0,1,0,0,0};
    static Matrix4 basis(H);
    Vector4 v(x0,s0,x1,s1);
    return basis*v;
}

void DrawHermiteCurve(HDC hdc, Vector& P0, Vector& T0, Vector& P1, Vector& T1, COLORREF color)
{
    //get a0, a1, a2, a3
    Vector4 xcoeff = GetHermiteCoeff(P0[0], T0[0], P1[0], T1[0]);
    Vector4 ycoeff = GetHermiteCoeff(P0[1], T0[1], P1[1], T1[1]);

    for (double t = 0; t <= 1; t += 0.001)
    {
        Vector4 vt;
        vt[3] = 1;
        // t t^2  t^3
        for (int i = 2; i >= 0; i--)
            vt[i] = vt[i + 1] * t;

        int x = round(DotProduct(xcoeff, vt));
        int y = round(DotProduct(ycoeff, vt));
        SetPixel(hdc, x, y, color);
    }
}


void DrawBezierCurve(HDC hdc,Vector& P0,Vector& P1,Vector& P2,Vector& P3,COLORREF c)
{
    /* t = 1/3
       x(0) = x1
       x`(0) = (x2-x1) / (1/3 - 0) = 3(x2 - x1)
       x(1) = x4
       x`(1) = (x4-x3) / (1 - 2/3) = 3(x4 - x3)


       y(0) = y1
       y`(0) = 3 (y2 - y1)
       y(1) = y4
       y`(1) = 3(y4 - y3)                                        */

    Vector T0(3*(P1[0]-P0[0]),3*(P1[1]-P0[1]));
    Vector T1(3*(P3[0]-P2[0]),3*(P3[1]-P2[1]));
    DrawHermiteCurve(hdc,P0,T0,P3,T1,c);
}



//draw square with Hermite curve
void FillingSquareWithHermite(HDC hdc, Vector& p1, double sideLen, COLORREF c)
{
    double left = p1[0];
    double right = p1[0] + sideLen;
    double top = p1[1];
    double bottom = p1[1] + sideLen;

    lineDDA(hdc, left, top, right, top, RGB(0, 0, 0));
    lineDDA(hdc, left, top, left, bottom, RGB(0, 0, 0));
    lineDDA(hdc, right, top, right, bottom, RGB(0, 0, 0));
    lineDDA(hdc, left, bottom, right, bottom, RGB(0, 0, 0));

    double controlLen = sideLen / 16;
    Vector T1(0, -controlLen);  // Top control point
    Vector T2(0, controlLen);   // Bottom control point
    for(left; left < right; left+=sideLen/500)
    {
        Vector P1(left, top);
        Vector P2(left, bottom);
        DrawHermiteCurve(hdc, P1, T1, P2, T2, c);
    }

}


void FillingRectangleWithBezier(HDC hdc, Vector& p1, double w, double h, COLORREF c)
{
    double left = p1[0];
    double right = p1[0] + w;
    double top = p1[1];
    double bottom = p1[1] + h;

    lineDDA(hdc, left, top, right, top, RGB(0, 0, 0));
    lineDDA(hdc, left, top, left, bottom, RGB(0, 0, 0));
    lineDDA(hdc, right, top, right, bottom, RGB(0, 0, 0));
    lineDDA(hdc, left, bottom, right, bottom, RGB(0, 0, 0));

    double stepSize = h / 500;
    for (double y = top; y < bottom; y += stepSize)
    {
        Vector T2(left + w / 100, y + h / 100);
        Vector T1(right - w / 100, y + h / 100);
        Vector P1(left, y);
        Vector P2(right, y);
        DrawBezierCurve(hdc, P1, T1, T2, P2, c);
    }
}




/*______________________________________ Filling algorithms ________________________________________*/


struct Entry
{
    int xmin,xmax;
};

struct EdgeRec
{
    double x;
    double minv;
    int ymax;
    bool operator<(EdgeRec r)
    {
        return x<r.x;
    }
};

typedef list<EdgeRec> EdgeList;
EdgeRec InitEdgeRec(POINT& v1,POINT& v2)
{
    if(v1.y>v2.y)swap(v1,v2);
    EdgeRec rec;
    rec.x=v1.x;
    rec.ymax=v2.y;
    rec.minv=(double)(v2.x-v1.x)/(v2.y-v1.y);
    return rec;
}
void InitEdgeTable(POINT *polygon,int n,EdgeList table[])
{
    POINT v1=polygon[n-1];
    for(int i=0; i<n; i++)
    {
        POINT v2=polygon[i];
        if(v1.y==v2.y)
        {
            v1=v2;
            continue;
        }
        EdgeRec rec=InitEdgeRec(v1, v2);
        table[v1.y].push_back(rec);
        v1=polygon[i];
    }
}
//nonConvex
void GeneralPolygonFill(HDC hdc,POINT *polygon,int n,COLORREF c)
{
    EdgeList *table=new EdgeList [MAXENTRIES];
    InitEdgeTable(polygon,n,table);
    int y=0;
    while(y<MAXENTRIES && table[y].size()==0)y++;

    if(y==MAXENTRIES)return;
    EdgeList ActiveList=table[y];
    while (ActiveList.size()>0)
    {
        ActiveList.sort();
        for(EdgeList::iterator it=ActiveList.begin(); it!=ActiveList.end(); it++)
        {
            int x1=(int)ceil(it->x);
            it++;
            int x2=(int)floor(it->x);
            for(int x=x1; x<=x2; x++)SetPixel(hdc,x,y,c);
        }
        y++;
        EdgeList::iterator it=ActiveList.begin();
        while(it!=ActiveList.end())
            if(y==it->ymax) it=ActiveList.erase(it);
            else it++;
        for(EdgeList::iterator it=ActiveList.begin(); it!=ActiveList.end(); it++)
            it->x+=it->minv;
        ActiveList.insert(ActiveList.end(),table[y].begin(),table[y].end());
    }
    delete[] table;
}


void InitEntries(Entry table[])
{
    for(int i=0; i<MAXENTRIES; i++)
    {
        table[i].xmin= INT_MAX;
        table[i].xmax= -INT_MAX;
    }
}

void ScanEdge(POINT v1,POINT v2,Entry table[])
{
    if(v1.y==v2.y)return;
    if(v1.y>v2.y)
        swap(v1,v2);
    double m=(double)(v2.x-v1.x)/(v2.y-v1.y);
    double x=v1.x;
    int y=v1.y;
    while(y<v2.y)
    {
        if(x<table[y].xmin)
            table[y].xmin=(int)ceil(x);
        if(x>table[y].xmax)
            table[y].xmax=(int)floor(x);
        y++;
        x+=m;
    }
}

void DrawScanLines(HDC hdc,Entry table[],COLORREF color)
{
    for(int y=0; y<MAXENTRIES; y++)
        if(table[y].xmin<table[y].xmax)
            for(int x=table[y].xmin; x<=table[y].xmax; x++)
                SetPixel(hdc,x,y,color);
}

//Convex
void ConvexFill(HDC hdc,POINT p[],int n,COLORREF color)
{
    Entry *table=new Entry[MAXENTRIES];
    InitEntries(table);
    POINT v1=p[n-1];
    for(int i=0; i<n; i++)
    {
        POINT v2=p[i];
        ScanEdge(v1,v2,table);
        v1=p[i];
    }
    DrawScanLines(hdc,table,color);
    delete table;
}


struct Vertex
{
    int x,y;
    Vertex(int x = 0,int y = 0):x(x),y(y)
    {
    }
};
//non recursive flood fill
void NRFloodFill(HDC hdc,int x,int y,COLORREF Cb,COLORREF Cf)
{
    stack<Vertex> S;
    Vertex v(x,y);
    S.push(v);
    while(!S.empty())
    {
        Vertex v1=S.top();
        S.pop();
        COLORREF c=GetPixel(hdc,v1.x,v1.y);
        if(c==Cb || c==Cf)continue;
        SetPixel(hdc,v1.x,v1.y,Cf);
        Vertex v2(v1.x, v1.y - 1);
        S.push(v2);
        Vertex v3(v1.x, v1.y + 1);
        S.push(v3);
        Vertex v4(v1.x + 1, v1.y);
        S.push(v4);
        Vertex v5(v1.x - 1, v1.y);
        S.push(v5);
    }
}
//recursive flood fill
void RFloodFill(HDC hdc,int x,int y,COLORREF Cb,COLORREF Cf )
{
    COLORREF C=GetPixel(hdc,x,y);
    if(C==Cb || C==Cf)return;
    SetPixel(hdc,x,y,Cf);
    RFloodFill(hdc,x+1,y,Cb,Cf );
    RFloodFill(hdc,x-1,y,Cb,Cf );
    RFloodFill(hdc,x,y+1,Cb,Cf);
    RFloodFill(hdc,x,y-1,Cb,Cf);
}

//Cardinal Spline curve.
void DrawCardinalSpline(HDC hdc, Vector P[], int n, double c, COLORREF C)
{
    double c1 = 1-c;
    Vector T0 = c1 * ((P[2][0] - P[0][0]), (P[2][1] - P[0][1]));
    for (int i = 0; i < n - 1; i++)
    {
        Vector T1 = c1 * ((P[i + 1][0] - P[i - 1][0]), (P[i + 1][1] - P[i - 1][1]));
        DrawHermiteCurve(hdc, P[i], T0, P[i + 1], T1, C);
        T0 = T1;
    }

}


/*______________________________________ Ellipse ________________________________________*/





void Draw4Points(HDC hdc, int xc, int yc, int x, int y, COLORREF color)
{
    SetPixel(hdc, xc + x, yc + y, color);  // First quadrant
    SetPixel(hdc, xc - x, yc + y, color);  // Second quadrant
    SetPixel(hdc, xc - x, yc - y, color);  // Third quadrant
    SetPixel(hdc, xc + x, yc - y, color);  // Fourth quadrant
}

void ellipseDirect(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
    int x, y;
    int a2 = a * a;
    int b2 = b * b;
    int twoA2 = 2 * a2;
    int twoB2 = 2 * b2;

    int xStart = xc - a;
    int xEnd = xc + a;

    for (x = xStart; x <= xEnd; x++)
    {
        y = yc + round(b * sqrt(1 - (double)(x - xc) * (x - xc) / a2));
        SetPixel(hdc, x, y, color);
        y = yc - round(b * sqrt(1 - (double)(x - xc) * (x - xc) / a2));
        SetPixel(hdc, x, y, color);
    }
}


void ellipsePolar(HDC hdc, int xc, int yc, int x0, int y0, COLORREF color)
{
    double dTheta = 1.0 / max(x0, y0);  // Step size for theta

    for (double theta = 0; theta <= 0.5*M_PI; theta += dTheta)
    {
        double x = x0 * cos(theta);
        double y = y0 * sin(theta);
        Draw4Points(hdc,xc, yc, (int)x, (int) y, color);
    }
}


void ellipseMidPoint(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
    int a2 = a * a;
    int b2 = b * b;
    int twoA2 = 2 * a2;
    int twoB2 = 2 * b2;
    int x = 0;
    int y = b;
    int dx = 0;
    int dy = twoA2 * y;

    // Region 1
    int p1 = b2 - a2 * b + (a2 / 4);

    while (dx < dy)
    {
        Draw4Points(hdc, xc, yc, x, y, color);

        if (p1 < 0)
        {
            x++;
            dx += twoB2;
            p1 += dx + b2;
        }
        else
        {
            x++;
            y--;
            dx += twoB2;
            dy -= twoA2;
            p1 += dx - dy + b2;
        }
    }

    // Region 2
    int p2 = b2 * (x + 0.5) * (x + 0.5) + a2 * (y - 1) * (y - 1) - a2 * b2;

    while (y >= 0)
    {
        Draw4Points(hdc, xc, yc, x, y, color);

        if (p2 > 0)
        {
            y--;
            dy -= twoA2;
            p2 += a2 - dy;
        }
        else
        {
            y--;
            x++;
            dx += twoB2;
            dy -= twoA2;
            p2 += dx - dy + a2;
        }
    }
}


/*______________________________________ Clipping with rectangle and Square ________________________________________*/


//Point
void PointClipping(HDC hdc,int x,int y,int xleft,int ytop,int xright,int ybottom,COLORREF color)
{
    if(x>=xleft && x<= xright && y>=ytop && y<=ybottom)
        SetPixel(hdc,x,y,color);
}


union OutCode
{
    unsigned All:4;
    struct
    {
        unsigned left:1,top:1,right:1,bottom:1;
    };
};
OutCode GetOutCode(double x,double y,int xleft,int ytop,int xright,int ybottom)
{
    OutCode out;
    out.All=0;
    if(x<xleft)out.left=1;
    else if(x>xright)out.right=1;
    if(y<ytop)out.top=1;
    else if(y>ybottom)out.bottom=1;
    return out;
}
void VIntersect(double xs,double ys,double xe,double ye,int x,double *xi,double *yi)
{
    *yi=ys+(x-xs)*(ye-ys)/(xe-xs);
    *xi=x;
}
void HIntersect(double xs,double ys,double xe,double ye,int y,double *xi,double *yi)
{
    *xi=xs+(y-ys)*(xe-xs)/(ye-ys);
    *yi=y;
}
//Line
void CohenSuth(HDC hdc,int xs,int ys,int xe,int ye,int xleft,int ytop,int xright,int ybottom)
{
    double x1=xs,y1=ys,x2=xe,y2=ye;
    OutCode out1=GetOutCode(x1,y1,xleft,ytop,xright,ybottom);
    OutCode out2=GetOutCode(x2,y2,xleft,ytop,xright,ybottom);
    while( (out1.All || out2.All) && !(out1.All & out2.All))
    {
        double xi,yi;
        if(out1.All)
        {
            if(out1.left)VIntersect(x1,y1,x2,y2,xleft,&xi,&yi);
            else if(out1.top)HIntersect(x1,y1,x2,y2,ytop,&xi,&yi);
            else if(out1.right)VIntersect(x1,y1,x2,y2,xright,&xi,&yi);
            else HIntersect(x1,y1,x2,y2,ybottom,&xi,&yi);
            x1=xi;
            y1=yi;
            out1=GetOutCode(x1,y1,xleft,ytop,xright,ybottom);
        }
        else
        {
            if(out2.left)VIntersect(x1,y1,x2,y2,xleft,&xi,&yi);
            else if(out2.top)HIntersect(x1,y1,x2,y2,ytop,&xi,&yi);
            else if(out2.right)VIntersect(x1,y1,x2,y2,xright,&xi,&yi);
            else HIntersect(x1,y1,x2,y2,ybottom,&xi,&yi);
            x2=xi;
            y2=yi;
            out2=GetOutCode(x2,y2,xleft,ytop,xright,ybottom);
        }
    }
    if(!out1.All && !out2.All)
    {
        MoveToEx(hdc,round(x1),round(y1),NULL);
        LineTo(hdc,round(x2),round(y2));
    }
}

//Polygon.
typedef vector<Vertex> VertexList;
typedef bool (*IsInFunc)(Vertex& v,int edge);
typedef Vertex (*IntersectFunc)(Vertex& v1,Vertex& v2,int edge);

VertexList ClipWithEdge(VertexList p,int edge,IsInFunc In,IntersectFunc Intersect)
{
    VertexList OutList;
    Vertex v1=p[p.size()-1];
    bool v1_in=In(v1,edge);
    for(int i=0; i<(int)p.size(); i++)
    {
        Vertex v2=p[i];
        bool v2_in=In(v2,edge);
        if(!v1_in && v2_in)
        {
            OutList.push_back(Intersect(v1,v2,edge));
            OutList.push_back(v2);
        }
        else if(v1_in && v2_in) OutList.push_back(v2);
        else if(v1_in) OutList.push_back(Intersect(v1,v2,edge));
        v1=v2;
        v1_in=v2_in;
    }
    return OutList;
}
bool InLeft(Vertex& v,int edge)
{
    return v.x>=edge;
}
bool InRight(Vertex& v,int edge)
{
    return v.x<=edge;
}
bool InTop(Vertex& v,int edge)
{
    return v.y>=edge;
}
bool InBottom(Vertex& v,int edge)
{
    return v.y<=edge;
}

Vertex VIntersect(Vertex& v1,Vertex& v2,int xedge)
{
    Vertex res;
    res.x=xedge;
    res.y=v1.y+(xedge-v1.x)*(v2.y-v1.y)/(v2.x-v1.x);
    return res;
}
Vertex HIntersect(Vertex& v1,Vertex& v2,int yedge)
{
    Vertex res;
    res.y=yedge;
    res.x=v1.x+(yedge-v1.y)*(v2.x-v1.x)/(v2.y-v1.y);
    return res;
}
void PolygonClip(HDC hdc,POINT *p,int n,int xleft,int ytop,int xright,int ybottom)
{
    VertexList vlist;
    for(int i=0; i<n; i++)vlist.push_back(Vertex(p[i].x,p[i].y));
    vlist=ClipWithEdge(vlist,xleft,InLeft,VIntersect);
    vlist=ClipWithEdge(vlist,ytop,InTop,HIntersect);
    vlist=ClipWithEdge(vlist,xright,InRight,VIntersect);
    vlist=ClipWithEdge(vlist,ybottom,InBottom,HIntersect);
    Vertex v1;
    if(vlist.size()>0)
    {
        v1=vlist[vlist.size()-1];
        for(int i=0; i<(int)vlist.size(); i++)
        {
            Vertex v2=vlist[i];
            MoveToEx(hdc,round(v1.x),round(v1.y),NULL);
            LineTo(hdc,round(v2.x),round(v2.y));
            v1=v2;
        }
    }

}


int x_left,x_right,y_bottom,y_top;
void drawRectangle(HDC hdc,int x1,int yy1,int d1,int d2,COLORREF c)
{
    x_left=x1;
    x_right=x1+d1;
    y_top=yy1;
    y_bottom=yy1+d2;
    lineDDA(hdc, x1, yy1, x_right,yy1, c);
    lineDDA(hdc, x1, yy1, x1,y_bottom, c);
    lineDDA(hdc, x1, y_bottom, x_right,y_bottom, c);
    lineDDA(hdc, x_right, yy1, x_right,y_bottom, c);

}
void drawSquare(HDC hdc,int x1,int yy1,int d,COLORREF c)
{
    x_left=x1;
    x_right=x1+d;
    y_top=yy1;
    y_bottom=yy1+d;
    lineDDA(hdc, x1, yy1, x_right,yy1, c);
    lineDDA(hdc, x1, yy1, x1,y_bottom, c);
    lineDDA(hdc, x1, y_bottom, x_right,y_bottom, c);
    lineDDA(hdc, x_right, yy1, x_right,y_bottom, c);
}


/*__________________________________________ Menu _______________________________________*/

HMENU list_;

void AddMenus(HWND hwnd)
{
    list_ = CreateMenu();
    HMENU color_list = CreateMenu();
    HMENU circle_list = CreateMenu();
    HMENU line_list = CreateMenu();
    HMENU Ellipse_list = CreateMenu();
    HMENU Filling_list = CreateMenu();
    HMENU Quarter_list = CreateMenu();
    HMENU clipping_list = CreateMenu();
    HMENU Clear_list = CreateMenu();
    HMENU Spiling_list = CreateMenu();

    /////////////////////////////////////
    AppendMenuW(line_list, MF_STRING, 1, L"DDA");
    AppendMenuW(line_list, MF_STRING, 2, L"MidPoint");
    AppendMenuW(line_list, MF_STRING, 3, L"Parametric");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)line_list, L"Lines");



    AppendMenuW(circle_list, MF_STRING, 4, L"Direct");
    AppendMenuW(circle_list, MF_STRING, 5, L"Polar");
    AppendMenuW(circle_list, MF_STRING, 6, L"Iterative Polar");
    AppendMenuW(circle_list, MF_STRING, 7, L"MidPoint");                    //circles
    AppendMenuW(circle_list, MF_STRING, 8, L" MidPoint Modification");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)circle_list, L"Circle");



    AppendMenuW(Quarter_list, MF_STRING, 9, L"1");
    AppendMenuW(Quarter_list, MF_STRING, 10, L"2");
    AppendMenuW(Quarter_list, MF_STRING, 11, L"3");
    AppendMenuW(Quarter_list, MF_STRING, 12, L"4");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)Quarter_list, L"Quarters");
    AppendMenuW(Filling_list, MF_STRING, 13, L"Filling quarter with lines");  // Fill the shape
    AppendMenuW(Filling_list, MF_STRING, 14, L"Filling quarter with circles");
    AppendMenuW(Filling_list, MF_STRING, 15, L"Filling Square with Hermit Curve");
    AppendMenuW(Filling_list, MF_STRING, 16, L"Filling Rectangle with Bezier Curve");


    AppendMenuW(Filling_list, MF_STRING, 17, L"flood fill - recursive"); //flood fill
    AppendMenuW(Filling_list, MF_STRING, 18, L"flood fill - non recursive"); //flood fill    // Filling
    AppendMenuW(Filling_list, MF_STRING, 19, L"Filling General Polygon");
    AppendMenuW(Filling_list, MF_STRING, 20, L"Filling Convex Polygon");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)Filling_list, L"Fillings");

    AppendMenuW(Spiling_list, MF_STRING, 21, L"Spiline");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)Spiling_list, L"Spiline");



    AppendMenuW(Ellipse_list, MF_STRING, 22, L"Direct");
    AppendMenuW(Ellipse_list, MF_STRING, 23, L"Polar");
    AppendMenuW(Ellipse_list, MF_STRING, 24, L"MidPoint");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)Ellipse_list, L"Ellipse");



    AppendMenuW(clipping_list, MF_STRING, 25, L"Clipping polygon with ractangle window");
    AppendMenuW(clipping_list, MF_STRING, 26, L"Clipping line with ractangle window");
    AppendMenuW(clipping_list, MF_STRING, 27, L"Clipping point with ractangle window");
    AppendMenuW(clipping_list, MF_STRING, 28, L"Clipping line with square window");
    AppendMenuW(clipping_list, MF_STRING, 29, L"Clipping point with square window");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)clipping_list, L"clipping");


    AppendMenuW(color_list, MF_STRING, 30, L"Red");
    AppendMenuW(color_list, MF_STRING, 31, L"Green");
    AppendMenuW(color_list, MF_STRING, 32, L"Blue");
    AppendMenuW(color_list, MF_STRING, 33, L"Black");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)color_list, L"Color");



    AppendMenuW(Clear_list, MF_STRING, 34, L"Clear");
    AppendMenuW(Clear_list, MF_STRING, 35, L"Save");
    AppendMenuW(Clear_list, MF_STRING, 36, L"Load");

    AppendMenuW(list_, MF_POPUP, (UINT_PTR)Clear_list, L"File");


    HMENU Mouse_cursor=CreateMenu();
    AppendMenuW(Mouse_cursor, MF_STRING, 37, L"Arrow");
    AppendMenuW(Mouse_cursor, MF_STRING, 38, L"Cross");
    AppendMenuW(Mouse_cursor, MF_STRING, 39, L"Hand");
    AppendMenuW(Mouse_cursor, MF_STRING, 40, L"Help");
    AppendMenuW(list_, MF_POPUP, (UINT_PTR)Mouse_cursor, L"Mouse cursor");


    SetMenu(hwnd, list_);
}



/*__________________________________________________ Window Procedure _____________________________________*/

/*  This function is called by the Windows function DispatchMessage()  */

COLORREF c = RGB(0,0,0);
int case_number = 0, quarter = 1;
LRESULT CALLBACK WindowProcedure (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    HDC hdc = GetDC(hwnd);
    static int x1, x2, x3 = 0,xf,xs,xe;
    static int y1, y2, y3, x4 = 0, y4 = 0,yf,ys,ye;
    static POINT P[6];
    static Vector v[6];
    static int cnt = 0,cnt2=0,cnt3=0;
    static Vector p1;  // First point of the shape
    static Vector p2;  // Second point of the shape
    static double r;


    switch (message)                  /* handle the messages */
    {

    case WM_SETCURSOR:
    SetCursor(LoadCursor(NULL,mouse));
    break;

    case WM_CREATE:
        AddMenus(hwnd);
        break;
    case WM_COMMAND:
        cnt = 0;
        switch (wParam)
        {
        case(1):
            case_number = 1;
            cout << "Draw Line using DDA\n";
            break;
        case(2):
            case_number = 2;
            cout << "Draw Line using MidPoint\n";
            break;
        case(3):
            case_number = 3;
            cout << "Draw Line using Parametric\n";
            break;
        case(4):
            case_number = 4;
            cout << "Draw Circle using Direct\n";
            break;
        case(5):
            case_number = 5;
            cout << "Draw Circle using Polar\n";
            break;
        case(6):
            case_number = 6;
            cout << "Draw Circle using Iterative Polar\n";
            break;
        case(7):
            case_number = 7;
            cout << "Draw Circle using MidPoint\n";
            break;
        case(8):
            case_number = 8;
            cout << "Draw Circle using Modification MidPoint\n";
            break;

        case  (9):
            quarter = 1;
            cout << "first quarter" << endl;
            break;
        case  (10):
            quarter = 2;
            cout << "second quarter" << endl;
            break;
        case  (11):

            quarter = 3;
            cout << "third quarter" << endl;
            break;

        case  (12):
            quarter = 4;
            cout << "fourth quarter" << endl;
            break;

        case  (13):
            case_number = 13;
            cout << "Filling the chosen quarter with lines"<<endl;
            break;
        case  (14):
            case_number = 14;
            cout << "Filling the chosen quarter with circles"<<endl;

            break;
        case  (15):
            case_number = 15;
            cout << "Filling Square with Hermit Curve"<<endl;;
            break;
        case  (16):
            case_number = 16;
            cout << "Filling Rectangle with Bezier Curve"<<endl;;
            break;
        case (17):
            case_number = 17;
            cout << "filling with flood fill recursive" << endl;
            break;

        case (18):
            case_number = 18;
            cout << "filling with flood fill non-recursive" << endl;
            break;

        case (19):
            case_number = 19;
            cout << "filling general polygon" << endl;
            break;

        case (20):
            case_number = 20;
            cout << "filling convex polygon" << endl;
            break;

        case (21):
            case_number = 21;
            cout << "Sipline" << endl;
            break;

        case (22):
            case_number = 22;
            cout << "Draw Ellipse using Direct Algorithm" << endl;
            break;

        case (23):
            case_number = 23;
            cout << "Draw Ellipse using polar Algorithm" << endl;
            break;

        case (24):
            case_number = 24;
            cout << "Draw Ellipse using midPoint Algorithm" << endl;
            break;

        case(25):
            case_number=25;
            cout << "Clipping algorithms using Rectangle as Clipping Window[Polygon]\n";
            break;

        case(26):
            case_number=26;
            cout << "Clipping algorithms using Rectangle as Clipping Window[Line]\n";
            break;

        case(27):
            case_number=27;
            cout << "Clipping algorithms using Rectangle as Clipping Window[Point]\n";
            break;

        case(28):
            case_number=28;
            cout << "Clipping algorithms using Square as Clipping Window[Line]\n";
            break;

        case(29):
            case_number=29;
            cout << "Clipping algorithms using Square as Clipping Window[Point]\n";
            break;

        case(30):
            c = RGB(255, 0, 0);
            cout<<"Red Color\n";
            break;

        case(31):
            c = RGB(0, 255, 0);
            cout << "Green color\n";
            break;

        case(32):
            c = RGB(0, 0, 255);
            cout << "Blue color\n";
            break;

        case(33):
            c = RGB(0, 0, 0);
            cout << "Black color\n";
            break;

        case(34):
            InvalidateRect(hwnd, NULL, TRUE);
            cout << "Clear is done.\n";
            break;

        case(35):
            save(hwnd);
            cout << "Save is done.\n";
            break;

        case(36):
            load(hwnd, hdc);
            cout << "load is done.\n";
            break;
        case(37):
            mouse=IDC_ARROW;
            break;
        case(38):
            mouse=IDC_CROSS;
            break;

        case(39):
            mouse=IDC_HAND;
            break;
        case(40):
            mouse=IDC_HELP;
            break;

        }

    case WM_LBUTTONDOWN:
        hdc = GetDC(hwnd);

        // two clicks
        if((case_number >= 1 && case_number <= 16) || (case_number >= 22 && case_number <= 24))
        {
            if(cnt == 0)cnt++;
            else if(cnt == 1)
            {
                x1 = LOWORD(lParam);
                y1 = HIWORD(lParam);
                cnt++;
            }
            else if(cnt==2)
            {
                x2 = LOWORD(lParam);
                y2 = HIWORD(lParam);
                if(case_number == 1)
                {
                    lineDDA(hdc, x1, y1, x2, y2, c);
                    ReleaseDC(hwnd, hdc);
                }

                else if(case_number==2)
                {
                    lineBresenham(hdc, x1, y1, x2, y2, c);
                }

                else if(case_number==3)
                {
                    ParametricLine(hdc, x1, y1, x2, y2, c);
                }
                else if(case_number==4)
                {
                    r = distance(x1, y1, x2, y2);
                    CircleDirect(hdc,x1,y1,r,c);
                }
                else if(case_number==5)
                {
                    r = distance(x1, y1, x2, y2);
                    cicrlepolar(hdc,x1,y1,r,c);
                }
                else if(case_number==6)
                {
                    r = distance(x1, y1, x2, y2);
                    CircleIterativePolar(hdc,x1,y1,r,c);
                }
                else if(case_number==7)
                {
                    r= distance(x1, y1, x2, y2);
                    Circlemidpoint(hdc,x1,y1,r,c);
                }
                else if(case_number==8)
                {
                    r= distance(x1, y1, x2, y2);
                    CircleModifiedMidpoint(hdc,x1,y1,r,c);
                }
                else if(case_number==13)
                {
                    r= distance(x1, y1, x2, y2);
                    fillCircleWithLines(hdc, x1, y1, r, c,quarter);
                }
                else if(case_number==14)
                {
                    r = distance(x1, y1, x2, y2);
                    CircleFillWithCircles(hdc, x1, y1, r, c,quarter);
                }
                else if(case_number == 15)
                {
                    p1 = Vector(x1,y1);
                    p2 = Vector(x2,y2);
                    // Calculate the width and height of the shape
                    double width = abs(p2[0] - p1[0]);
                    double height = abs(p2[1] - p1[1]);
                    double sideLen = max(width, height);
                    FillingSquareWithHermite(hdc, p1, sideLen, c);
                }
                else if(case_number == 16)
                {
                    p1 = Vector(x1,y1);
                    p2 = Vector(x2,y2);
                    // Calculate the width and height of the shape
                    double width = abs(p2[0] - p1[0]);
                    double height = abs(p2[1] - p1[1]);
                    FillingRectangleWithBezier(hdc, p1, width, height, c);

                }
                else if(case_number == 22)
                {
                    int  majorAxis = distance(x1, y1, x2, y2);
                    int  minorAxis = distance(x1, y1, x1 + majorAxis / 2, y1);
                    ellipseDirect(hdc, x1, y1, majorAxis,minorAxis, c);
                }
                else if(case_number == 23)
                {
                    int  majorAxis = distance(x1, y1, x2, y2);
                    int  minorAxis = distance(x1, y1, x1 + majorAxis / 2, y1);
                    ellipsePolar(hdc, x1, y1, majorAxis,minorAxis, c);
                }
                else if(case_number == 24)
                {
                    int  majorAxis = distance(x1, y1, x2, y2);
                    int  minorAxis = distance(x1, y1, x1 + majorAxis / 2, y1);
                    ellipseMidPoint(hdc, x1, y1, majorAxis,minorAxis, c);
                }
                cnt=1;
            }
        }

        else if(case_number >= 17 && case_number <= 20)
        {
            if (cnt <= 8)
            {
                P[cnt-1].x = LOWORD(lParam);
                P[cnt-1].y = HIWORD(lParam);
                cnt++;
                if (cnt == 7)
                {
                    Polygon(hdc, P, 6);
                }
                if(cnt == 8)
                {
                    if(case_number==17)
                    {
                        xf=LOWORD(lParam);
                        yf = HIWORD(lParam);
                        RFloodFill(hdc,xf,yf,c,c);
                        ReleaseDC(hwnd, hdc);
                    }
                    if(case_number==18)
                    {
                        xf=LOWORD(lParam);
                        yf = HIWORD(lParam);
                        NRFloodFill(hdc,xf,yf,c,c);
                        ReleaseDC(hwnd, hdc);
                    }
                    if(case_number==19)
                    {
                        GeneralPolygonFill(hdc,P,6,c);
                        ReleaseDC(hwnd, hdc);
                    }
                    if(case_number==20)
                    {
                        ConvexFill(hdc, P, 6, c);
                        ReleaseDC(hwnd, hdc);
                    }
                    cnt=1;
                }
            }
        }

        else if(case_number==21){

          if (cnt <=7){

            int x= LOWORD(lParam);
            int y= HIWORD(lParam);
            v[cnt-1]=Vector(x,y);
            cnt++;
            if (cnt == 7)
            {
                DrawCardinalSpline(hdc,v,6,0.5,c);
                cnt=1;
                ReleaseDC(hwnd, hdc);
            }
        }
    }



        else if(case_number>=25 && case_number<=29)
        {
            if(cnt==0)cnt++;
            else if(cnt==1 && (case_number==29 ||case_number==27)) //point ->square,Rectangle
            {
                xs = LOWORD(lParam);
                ys = HIWORD(lParam);
                PointClipping(hdc,xs,ys,x_left,y_top,x_right,y_bottom,c);
                ReleaseDC(hwnd,hdc);
            }
            else if(case_number==28 || case_number==26) //line ->square,Rectangle
            {
                if(cnt==1)
                {
                    xs = LOWORD(lParam);
                    ys = HIWORD(lParam);
                    cnt++;
                }
                else if(cnt==2)
                {
                    xe = LOWORD(lParam);
                    ye = HIWORD(lParam);
                    CohenSuth(hdc,xs,ys,xe,ye,x_left,y_top,x_right,y_bottom);
                    ReleaseDC(hwnd,hdc);
                    cnt=1;
                }
            }
            else if(cnt==1 && case_number==25) //polygon ->Rectangle
            {
                static POINT points[5];
                if (cnt3 < 6)
                {
                    points[cnt3].x = LOWORD(lParam);
                    points[cnt3].y = HIWORD(lParam);
                    cnt3++;
                    if (cnt3 == 5)
                    {
                        PolygonClip(hdc,points,5,x_left,y_top,x_right,y_bottom);
                        ReleaseDC(hwnd, hdc);
                        cnt3=0;
                        cnt=1;
                    }
                }

            }
        }

        break;
    case WM_RBUTTONDOWN:
        hdc = GetDC(hwnd);
        if(case_number>=25 && case_number<=29)
        {
            if(cnt2==0)
            {
                x1 = LOWORD(lParam);
                y1 = HIWORD(lParam);
                cnt2++;
            }
            else if(cnt2==1)
            {
                x2 = LOWORD(lParam);
                y2 = HIWORD(lParam);
                if(case_number>=28 && case_number<=29)
                {
                    int d=distance(x1,y1,x2,y2);
                    drawSquare(hdc,x1,y1,d,c);
                    ReleaseDC(hwnd,hdc);
                    cnt2=0;
                }
                else cnt2++;
            }
            else if(case_number>=25 && case_number<=27 && cnt2==2) //rectangle
            {
                x3 = LOWORD(lParam);
                y3 = HIWORD(lParam);
                int d1=distance(x1,y1,x2,y2);
                int d2=distance(x2,y2,x3,y3);
                drawRectangle(hdc,x1,y1,d1,d2,c);
                ReleaseDC(hwnd,hdc);
                cnt2=0;
            }
        }
        break;

    case WM_DESTROY:
        PostQuitMessage(0);       /* send a WM_QUIT to the message queue */
        break;
    default:                      /* for messages that we don't deal with */
        return DefWindowProc (hwnd, message, wParam, lParam);
    }

    return 0;
}





/*________________________________ Save and load _______________________________________*/



bool HDCToFile(const char* FilePath, HDC Context, RECT Area, uint16_t BitsPerPixel = 24)
{
    uint32_t Width = Area.right - Area.left;
    uint32_t Height = Area.bottom - Area.top;

    BITMAPINFO Info;
    BITMAPFILEHEADER Header;
    memset(&Info, 0, sizeof(Info));
    memset(&Header, 0, sizeof(Header)); // bnsfrhom

    Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER); //set el size
    Info.bmiHeader.biWidth = Width; // set el width
    Info.bmiHeader.biHeight = Height; //set el height
    Info.bmiHeader.biPlanes = 1; //1 for non-multilayered images.
    Info.bmiHeader.biBitCount = BitsPerPixel; //number of bits per pixel in the image.
    Info.bmiHeader.biCompression = BI_RGB; //uncompressed
    Info.bmiHeader.biSizeImage = Width * Height * (BitsPerPixel > 24 ? 4 : 3); //24=>> 8+8+8 , r+g+b
    Header.bfType = 0x4D42; //0x4D42 represents the ASCII characters 'BM'
    Header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
    char* Pixels = NULL;
    HDC MemDC = CreateCompatibleDC(Context);
    HBITMAP Section = CreateDIBSection(Context, &Info, DIB_RGB_COLORS, (void**)&Pixels, 0, 0);
    DeleteObject(SelectObject(MemDC, Section));
    BitBlt(MemDC, 0, 0, Width, Height, Context, Area.left, Area.top, SRCCOPY);
    DeleteDC(MemDC);
    std::fstream hFile(FilePath, std::ios::out | std::ios::binary);
    if (hFile.is_open())
    {
        hFile.write((char*)&Header, sizeof(Header));
        hFile.write((char*)&Info.bmiHeader, sizeof(Info.bmiHeader));
        hFile.write(Pixels, (((BitsPerPixel * Width + 31) & ~31) / 8) * Height);
        hFile.close();
        DeleteObject(Section);
        return true; //saved
    }
    DeleteObject(Section);
    return false; //failed
}
void load(HWND hWnd, HDC &hdc)
{
    string fileName = "picture.bmp";
    if (fileName == "")
        return ;
    HBITMAP hBitmap; // handle to the actual bitmap image
    hBitmap = (HBITMAP)::LoadImage(NULL, fileName.c_str(), IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
    HDC hLocalDC;
    hLocalDC = CreateCompatibleDC(hdc);
    BITMAP qBitmap; // save retrieved data && holds information about the bitmap
    int iReturn = GetObject(reinterpret_cast<HGDIOBJ>(hBitmap), sizeof(BITMAP),reinterpret_cast<LPVOID>(&qBitmap)); // retrieve el data w nsavhaa fl qbitmap
    HBITMAP hOldBmp = (HBITMAP)SelectObject(hLocalDC, hBitmap);
    BOOL qRetBlit = BitBlt(hdc, 0, 0, qBitmap.bmWidth, qBitmap.bmHeight,hLocalDC, 0, 0, SRCCOPY);
    SelectObject (hLocalDC, hOldBmp);
    DeleteDC(hLocalDC);
    DeleteObject(hBitmap);
}
void save(HWND &hWnd)
{
    HDC hdc = GetDC(hWnd);
    string fileName = "picture.bmp";
    if (fileName == "")
        return ;
    int windowWidth ;
    int windowHeight;
    RECT rect; //The Rect class holds the four borders of a rectangle:
    // left, top, right, and bottom.
    if(GetWindowRect(hWnd, &rect))
    {
        windowWidth = rect.right - rect.left;
        windowHeight = rect.bottom - rect.top;
    }
    RECT rect1 = {0, 0, windowWidth, windowHeight};
    HDCToFile(fileName.c_str(),hdc,rect1);
}


