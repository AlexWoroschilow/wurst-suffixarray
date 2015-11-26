#include <CoreFoundation/CoreFoundation.h>
#include <CoreServices/CoreServices.h>
#include <QuickLook/QuickLook.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <fcntl.h>
#include "e_malloc.h"
#include <string.h>
#include <math.h>
#include <GLUT/glut.h>
#include <OPENGL/glu.h> 
#include <OPENGL/glext.h> 

//#include "main.h"
#include "specktacular.h"
#include "coord.h"
#include "pdbin_i.h"

#define bool int
#define false 0
#define true 1

#define MAXGRID 90
float sphi=90.0, stheta=45.0;
float sdepth = 10;
float zNear=0.0001, zFar=500000.0;
float aspect = 4/3;//kWindowWidth/kWindowHeight;
long xsize, ysize;
int downX, downY;
bool leftButton = false, middleButton = false;
struct speckstate *state;
GLuint dl;
enum{BACKBONE_MODE, TRACE_MODE, RIBBON_MODE, RIBBON3D_MODE};
const int MAX_DRAW_MODES = 4;
enum{PPL2_SHADER, PHONG_SHADER, GOOCH_SHADER, LATTICE_SHADER, TEST_SHADER};
const int MAX_SHADE_MODES = 5;
int drawmode = RIBBON_MODE;
int shademode = PHONG_SHADER;


/*
 *  from main.c
 *  gltest
 *
 *  Created by Thomas Margraf on 21/6/2007.
 *  Copyright 2007 Thomas Margraf. All rights reserved.
 *
 */

/* TODO:
 - rollercoster mode
 - 3D-ribbon
 - rocking mode
 - webcam texturing
 - antialiasing/jitter via accumulation buffer if !mousedown
 - surfaces (van-der-waal's surf, SAS)
 - transparency
 - high voltage H-bonds - maybe as a particle system in vertex shader.
 - detect chain breaks (check bond lengths)
 - wiimote control
 - support RNA/DNA
 - support selection/highlighting of multiple regions
 */

/***************************************************************************
 readShaderSrc - read in an ASCII file with GLSL sources
 ****************************************************************************/
char *readShaderSrc(char *fn) {
    FILE *fp;
    char *content = NULL;
	
    int count=0;
	
    if (fn != NULL) {
        fp = fopen(fn,"rt");
        if (fp != NULL) {
            fseek(fp, 0, SEEK_END);
            count = ftell(fp);
            rewind(fp);
			
            if (count > 0) {
                content = (char *)E_MALLOC(sizeof(char) * (count+1));
                count = fread(content,sizeof(char),count,fp);
                content[count] = '\0';
            }
            fclose(fp);
        }
    }
    return content;
}


/***************************************************************************
 cross - Calculate the cross product and return it
 ****************************************************************************/
static void cross (float dst[3], float srcA[3], float srcB[3])
{
    dst[0] = srcA[1]*srcB[2] - srcA[2]*srcB[1];
    dst[1] = srcA[2]*srcB[0] - srcA[0]*srcB[2];
    dst[2] = srcA[0]*srcB[1] - srcA[1]*srcB[0];
}


/***************************************************************************
 normalize - Normalize the input vector
 ****************************************************************************/
static void normalize (float vec[3])
{
    const float squaredLen = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
    const float invLen = 1.f / (float) sqrt (squaredLen);
    
    vec[0] *= invLen;
    vec[1] *= invLen;
    vec[2] *= invLen;
}


/***************************************************************************
 scale - Scale the given vector
 ****************************************************************************/
static void scale (float v[3], float s)
{
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
}

/***************************************************************************
 dist
 ****************************************************************************/
static float dist(GLfloat* ctrlpts, int i, int j){
    float dx, dy, dz;
    dx = ctrlpts[i]-ctrlpts[j];
    dy = ctrlpts[i+1]-ctrlpts[j+1];
    dz = ctrlpts[i+2]-ctrlpts[j+2];
    dx *= dx;
    dy *= dy;
    dz *= dz;
    return(sqrt(dx+dy+dz));
}


/***************************************************************************
 DrawProt
 ****************************************************************************/
static void DrawProt(int snum){
    size_t i;
    float m[16];
    float *xaxis = &m[0],
	*up = &m[4],
	*at = &m[8];
    float comx = 0.0,
	comy = 0.0,
	comz = 0.0;
    struct coord *c = state->coords[snum];
    for (i=0; i<c->size; i++){
        comx += c->rp_ca[i].x;
        comy += c->rp_ca[i].y;
        comz += c->rp_ca[i].z;
    }
    comx /= c->size;
    comy /= c->size;
    comz /= c->size;
	
    
    GLUquadricObj *qobj;
    
    qobj = gluNewQuadric();
    gluQuadricDrawStyle(qobj, GLU_FILL);
    gluQuadricNormals(qobj, GLU_SMOOTH);
    
    gluQuadricTexture(qobj, GL_TRUE);
    
    glColor3f(0.0, 0.0, 1.0);
    for (i=0; i<c->size; i++){ 
        glPushMatrix();
        glTranslatef(c->rp_ca[i].x - comx, c->rp_ca[i].y - comy, c->rp_ca[i].z - comz);
        gluSphere(qobj, 0.5, 8, 8);
        glPopMatrix();
    }
    
    
    glColor4f(0.8, 0.8, 0.8, 1.0);
    
    if(snum%3 == 1){
        glColor3f(1.0, 0.0, 0.0);
    }
    else if(snum%3 == 2){
        glColor3f(0.0, 1.0, 0.0);
        
    }
    else{
        glColor3f(0.0, 1.0, 0.0);
    }
    
    for (i=1; i<c->size; i++){
        float dx, dy, dz;
        float mx, my, mz;
        float len;
		
        mx = c->rp_ca[i].x - comx;
        my = c->rp_ca[i].y - comy;
        mz = c->rp_ca[i].z - comz;
        glPushMatrix();
        glTranslatef(mx, my, mz);
        
        dx = c->rp_ca[i].x - c->rp_ca[i-1].x;
        dy = c->rp_ca[i].y - c->rp_ca[i-1].y;
        dz = c->rp_ca[i].z - c->rp_ca[i-1].z;
        len = sqrt(dx*dx + dy*dy + dz*dz);
        if(len < 4.9){
            at[0]= dx;
            at[1]= dy;
            at[2]= dz;
            normalize (at);
            // Make a useable copy of the current up vector.
            up[0] = 0.0; up[1] = 1.0; up[2] = 0.0;
            
            // Cross product of the new look at vector and the current
            //   up vector will produce a vector which is the new
            //   positive X axis of our transformed object.
            cross (xaxis, at, up);
            normalize (xaxis);
            
            // Calculate the new up vector, which will be the
            //   positive Y axis of our transformed object. Note
            //   that it will lie in the same plane as the new
            //   look at vector and the old up vector.
            cross (up, xaxis, at);
            
            // Account for the fact that the geometry will be defined to
            //   point along the negative Z axis.
            scale (at, -1.f);
            
            // Fill out the rest of the 4x4 matrix
            m[3] = 0.f;     // xaxis is m[0..2]
            m[7] = 0.f;     // up is m[4..6]
            m[11] = 0.f;    // -at is m[8..10]
            m[12] = 0.0; m[13] = 0.0; m[14] = 0.0;
            m[15] = 1.f;
            
            // Multiply onto current matrix stack.
            glPushMatrix();
            glMultMatrixf(m);
            gluCylinder(qobj, 0.5, 0.5, len, 8, 1);
            glPopMatrix();
        }
        glPopMatrix();
    }
    
}

void nurbsError(GLenum errorCode){
    const GLubyte *estring;
    estring = gluErrorString(errorCode);
    fprintf(stderr, "NURBS error: %s\n", estring);
}

/***************************************************************************
 DrawRibbon
 ****************************************************************************/
static void DrawRibbon(int snum){
    int i = 0;
    float comx = 0.0,
	comy = 0.0,
	comz = 0.0;
    struct coord *c = state->coords[snum];
	
    for (i=0; i<c->size; i++){
        comx += c->rp_ca[i].x;
        comy += c->rp_ca[i].y;
        comz += c->rp_ca[i].z;
    }
    comx /= c->size;
    comy /= c->size;
    comz /= c->size;
    
    GLUnurbsObj *nurbs = gluNewNurbsRenderer();
    
    GLfloat* ctrlpts = E_MALLOC(6*sizeof(GLfloat)*c->size);
    for(i=0; i < c->size*6; i+=6){
        
		ctrlpts[i]   = c->rp_n[i/6].x - comx;
		ctrlpts[i+1] = c->rp_n[i/6].y - comy;
		ctrlpts[i+2] = c->rp_n[i/6].z - comz;
        
		ctrlpts[i+3] = c->rp_c[i/6].x - comx;
		ctrlpts[i+4] = c->rp_c[i/6].y - comy;
		ctrlpts[i+5] = c->rp_c[i/6].z - comz;
		
    }
    
    GLfloat* uknots = E_MALLOC(sizeof(GLfloat)*(2*c->size+4));
    
    for(i=0; i < 2*c->size+4; i++){
        if(i>c->size*2+3){
            uknots[i]=uknots[c->size*2+3];
        }
        else{
            uknots[i]=i/1.0; 
        }
    }
	
	GLfloat vknots[8] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
	
	glColor3f(0.0, 0.0, 1.0);
	glLineWidth(2);
	glEnable(GL_AUTO_NORMAL);
	glShadeModel(GL_SMOOTH);
	
	glPushMatrix();
	
	gluNurbsCallback(nurbs, GLU_ERROR, (GLvoid(*)())nurbsError);
	gluNurbsProperty(nurbs, GLU_V_STEP, 4);
	gluNurbsProperty(nurbs, GLU_U_STEP, 4);
	gluNurbsProperty(nurbs, GLU_CULLING, GL_TRUE);
	gluNurbsProperty(nurbs, GLU_SAMPLING_METHOD, GLU_DOMAIN_DISTANCE);
	
	glShadeModel(GL_SMOOTH);
	
	gluBeginSurface(nurbs);
	
	
	if(snum%3 == 1){
		glColor3f(1.0, 0.0, 0.0);
	}
	else if(snum%3 == 2){
		glColor3f(0.0, 1.0, 0.0);
		
	}
	else{
		glColor3f(0.0, 1.0, 0.0);
	}
	
	gluNurbsSurface(nurbs,              //context object
					c->size*2+4,          //number of u-knots
					uknots,                //u-knot pointer
					8,                  //number of v-knots
					vknots,               //v-knot ptr
					3,                  //width of u-control points (u-stride)
					3,                  //width of v-control points (v-stride)
					ctrlpts,  //control point pointer
					4,                  // u-order of the curve (degree+1)
					4,                  // v-order of the curve (degree+1)
					GL_MAP2_TEXTURE_COORD_2    //type
					);
	
	gluNurbsSurface(nurbs,              //context object
					c->size*2+4,          //number of u-knots
					uknots,                //u-knot pointer
					8,                  //number of v-knots
					vknots,               //v-knot ptr
					3,                  //width of u-control points (u-stride)
					3,                  //width of v-control points (v-stride)
					ctrlpts,  //control point pointer
					4,                  // u-order of the curve (degree+1)
					4,                  // v-order of the curve (degree+1)
					GL_MAP2_VERTEX_3    //type
					);
	
	gluEndSurface(nurbs);
	
	glPopMatrix();
}

/***************************************************************************
 Draw3dRibbon
 ****************************************************************************/
static void Draw3dRibbon(int snum){
    int i = 0;
    float comx = 0.0,
	comy = 0.0,
	comz = 0.0;
    float *cp = E_MALLOC(3*sizeof(float));
    float *inA = E_MALLOC(3*sizeof(float));
    float *inB = E_MALLOC(3*sizeof(float));
    float *inC = E_MALLOC(3*sizeof(float));
    float *out = E_MALLOC(3*sizeof(float));
	
    struct coord *c = state->coords[snum];
	
    for (i=0; i<c->size; i++){
        comx += c->rp_ca[i].x;
        comy += c->rp_ca[i].y;
        comz += c->rp_ca[i].z;
    }
    comx /= c->size;
    comy /= c->size;
    comz /= c->size;
    
    GLUnurbsObj *nurbs = gluNewNurbsRenderer();
    
    GLfloat* ctrlpts = E_MALLOC(12*sizeof(GLfloat)*(c->size-1));
    
    GLfloat* tracepts = E_MALLOC(3*sizeof(GLfloat)*c->size-1);
    
    for(i=0; i < c->size-1; i++){
        if((i) %2){
            ctrlpts[i*12]   = c->rp_c[i].x - comx;
            ctrlpts[i*12+1] = c->rp_c[i].y - comy;
            ctrlpts[i*12+2] = c->rp_c[i].z - comz;
        }
        else{
            ctrlpts[i*12]   = c->rp_o[i].x - comx;
            ctrlpts[i*12+1] = c->rp_o[i].y - comy;
            ctrlpts[i*12+2] = c->rp_o[i].z - comz;
        }
        
        inA[0] = c->rp_c[i].x - c->rp_ca[i].x;
        inA[1] = c->rp_c[i].y - c->rp_ca[i].y;
        inA[2] = c->rp_c[i].z - c->rp_ca[i].z;
        inB[0] = c->rp_c[i].x - c->rp_o[i].x;
        inB[1] = c->rp_c[i].y - c->rp_o[i].y;
        inB[2] = c->rp_c[i].z - c->rp_o[i].z;
        
        cross(cp, inA, inB);
        normalize(cp);
        
        inC[0] = ((c->rp_c[i].x + c->rp_ca[i].x) / 2) - comx;
        inC[1] = ((c->rp_c[i].y + c->rp_ca[i].y) / 2) - comy;
        inC[2] = ((c->rp_c[i].z + c->rp_ca[i].z) / 2) - comz;
        
        //cross(out, cp, inC);
        //normalize(out);
		
        ctrlpts[i*12+3] = inC[0] + cp[0];
        ctrlpts[i*12+4] = inC[1] + cp[1];
        ctrlpts[i*12+5] = inC[2] + cp[2];
        
        if((i) %2){
            ctrlpts[i*12+6] = c->rp_o[i].x - comx;
            ctrlpts[i*12+7] = c->rp_o[i].y - comy;
            ctrlpts[i*12+8] = c->rp_o[i].z - comz;
        }
        else{
            ctrlpts[i*12+6]   = c->rp_c[i].x - comx;
            ctrlpts[i*12+7] = c->rp_c[i].y - comy;
            ctrlpts[i*12+8] = c->rp_c[i].z - comz;
        }
        
        
		ctrlpts[i*12+9]  = inC[0] - cp[0];
		ctrlpts[i*12+10] = inC[1] - cp[1];
		ctrlpts[i*12+11] = inC[2] - cp[2];
		
		tracepts[i] = inC[0];
		tracepts[i+1] = inC[1];
		tracepts[i+2] = inC[2];
		
		//                printf("[%f, %f, %f], [%f, %f, %f] \n", 
		//                       cp[0], cp[1], cp[2],
		//                       inC[0], inC[1], inC[2]);
        
		//        printf("[%f, %f, %f], [%f, %f, %f], [%f, %f, %f], [%f, %f, %f] \n",
		//               ctrlpts[i], ctrlpts[i+1] ,ctrlpts[i+2], ctrlpts[i+3], ctrlpts[i+4],
		//               ctrlpts[i+5], ctrlpts[i+6], ctrlpts[i+7], ctrlpts[i+8], ctrlpts[i+9], 
		//               ctrlpts[i+10], ctrlpts[i+11]);
        
    }
	
    free(cp);
    free(inA);
    free(inB);
    free(inC);
    free(out);
    
	GLfloat* uknots = E_MALLOC(sizeof(GLfloat)*(c->size+3));
	
	for(i=0; i < c->size+3; i++){
		if(i<c->size+1 && i>2){
			uknots[i] = uknots[i-1]+1.0;
		}
		else if(i<3){
			uknots[i] = 0.0;
		}
		else{
			uknots[i] = uknots[i-1]; 
		}
	}
	
	
	GLfloat vknots[8] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
	
	GLfloat* knots = E_MALLOC(sizeof(GLfloat)*(c->size+4));
	
	for(i=0; i < c->size+4; i++){
		
		if(i>c->size+2){
			knots[i]=knots[c->size+2];
		}
		else{
			knots[i]=i/1.0-i%1; 
		}
	}
	
	
	
	glColor3f(0.0, 0.0, 1.0);
	glEnable(GL_AUTO_NORMAL);
	glShadeModel(GL_SMOOTH);
	
	glPushMatrix();
	
	gluNurbsCallback(nurbs, GLU_ERROR, (GLvoid(*)())nurbsError);
	gluNurbsProperty(nurbs, GLU_V_STEP, 4);
	gluNurbsProperty(nurbs, GLU_U_STEP, 4);
	gluNurbsProperty(nurbs, GLU_CULLING, GL_TRUE);
	gluNurbsProperty(nurbs, GLU_SAMPLING_METHOD, GLU_DOMAIN_DISTANCE);
	
	glShadeModel(GL_SMOOTH);
	
	gluBeginSurface(nurbs);
	
	gluNurbsSurface(nurbs,              //context object
					c->size+3,          //number of u-knots
					uknots,             //u-knot pointer
					8,                  //number of v-knots
					vknots,             //v-knot ptr
					3,                  //width of u-control points (u-stride)
					3,                  //width of v-control points (v-stride)
					ctrlpts,            //control point pointer
					4,                  // u-order of the curve (degree+1)
					4,                  // v-order of the curve (degree+1)
					GL_MAP2_TEXTURE_COORD_2    //type
					);
	
	gluNurbsSurface(nurbs,              //context object
					c->size+3,          //number of u-knots
					uknots,             //u-knot pointer
					8,                  //number of v-knots
					vknots,             //v-knot ptr
					3,                  //width of u-control points (u-stride)
					3,                  //width of v-control points (v-stride)
					ctrlpts,            //control point pointer
					4,                  // u-order of the curve (degree+1)
					4,                  // v-order of the curve (degree+1)
					GL_MAP2_NORMAL      //type
					);
	
	
	gluNurbsSurface(nurbs,              //context object
					c->size+3,          //number of u-knots
					uknots,             //u-knot pointer
					8,                  //number of v-knots
					vknots,             //v-knot ptr
					3,                  //width of u-control points (u-stride)
					3,                  //width of v-control points (v-stride)
					ctrlpts,            //control point pointer
					4,                  // u-order of the curve (degree+1)
					4,                  // v-order of the curve (degree+1)
					GL_MAP2_VERTEX_3    //type
					);
	
	
	gluEndSurface(nurbs);
	/*
	 gluBeginCurve(nurbs);
	 gluNurbsCurve(nurbs,              //context object
	 c->size+4,          //number of knots
	 knots,              //knot pointer
	 3,                  //width of control points
	 tracepts,  //control point pointer
	 4,                  //order of the curve (degree+1)
	 GL_MAP1_VERTEX_3    //type
	 );
	 
	 gluEndCurve(nurbs);
	 */
	glPopMatrix();
}


/***************************************************************************
 DrawTrace
 ****************************************************************************/
static void DrawTrace(int snum){
    int i = 0;
    float comx = 0.0,
	comy = 0.0,
	comz = 0.0;
    
    struct coord *c = state->coords[snum];
    
    for (i=0; i<c->size; i++){
        comx += c->rp_ca[i].x;
        comy += c->rp_ca[i].y;
        comz += c->rp_ca[i].z;
    }
    comx /= c->size;
    comy /= c->size;
    comz /= c->size;
    
    GLUnurbsObj *nurbs = gluNewNurbsRenderer();
    
    GLfloat* ctrlpts = E_MALLOC(12*sizeof(GLfloat)*c->size);
    for(i=0; i < c->size*12; i+=12){
        
        ctrlpts[i]   = c->rp_n[i/12].x - comx;
        ctrlpts[i+1] = c->rp_n[i/12].y - comy;
        ctrlpts[i+2] = c->rp_n[i/12].z - comz;
        
        ctrlpts[i+3] = 1.0;
        
        ctrlpts[i+4] = c->rp_ca[i/12].x - comx;
        ctrlpts[i+5] = c->rp_ca[i/12].y - comy;
        ctrlpts[i+6] = c->rp_ca[i/12].z - comz;
		
        ctrlpts[i+7] = 1.0;
        
        ctrlpts[i+8] = c->rp_c[i/12].x - comx;
        ctrlpts[i+9] = c->rp_c[i/12].y - comy;
        ctrlpts[i+10] = c->rp_c[i/12].z - comz;
		
        ctrlpts[i+11] = 1.0;
		
    }
    
    GLfloat* knots = E_MALLOC(sizeof(GLfloat)*(c->size*3+4));
    
    for(i=0; i < c->size*3+4; i++){
		
        if(i>c->size*3+2){
            knots[i]=knots[c->size*3+2];
        }
        else{
            knots[i]=i/1.0-i%1; 
        }
    }
    
    glColor3f(0.0, 0.0, 1.0);
    if(snum%3 == 1){
        glColor3f(1.0, 0.0, 0.0);
    }
    else if(snum%3 == 2){
        glColor3f(0.0, 1.0, 0.0);
        
    }
    else{
        glColor3f(0.0, 1.0, 0.0);
    }
    
    
    glLineWidth(2);
    
    glPushMatrix();
    
    gluNurbsCallback(nurbs, GLU_ERROR, (GLvoid(*)())nurbsError);
    
    gluBeginCurve(nurbs);
	
    gluNurbsCurve(nurbs,              //context object
                  c->size*3+4,          //number of knots
                  knots,              //knot pointer
                  4,                  //width of control points
                  ctrlpts,  //control point pointer
                  4,                  //order of the curve (degree+1)
                  GL_MAP1_VERTEX_4    //type
				  );
	
    gluEndCurve(nurbs);
    
    glPopMatrix();
}

/***************************************************************************
 InitGL - initialize the OpenGL system
 ****************************************************************************/
GLvoid InitGL(GLvoid)           // All Setup For OpenGL Goes Here 
{ 
    int i;
    float r, g, b = 1.0;
    GLfloat mat_specular[] = {r, g, b, 1.0};
    GLfloat mat_shininess[] = {50.0};
    //GLfloat light_position[] = {0.0, 0.0, 400.0, 0.0};
    //GLfloat white_light[] = {1.0, 1.0, 1.0, 1.0};
    //GLfloat lmodel_ambient[] = {0.3, 0.3, 0.3, 1.0};
    
    glShadeModel(GL_SMOOTH);        // Enables Smooth Shading 
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);    // Black Background
    glClearDepth(1.0f);          // Depth Buffer Setup 
    glEnable(GL_DEPTH_TEST);        // Enables Depth Testing 
    glDepthFunc(GL_LEQUAL);         // The Type Of Depth Test To Do 
	// Really Nice Perspective Calculations 
    //glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST); 
    
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glMaterialfv(GL_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_BACK, GL_SHININESS, mat_shininess);    
    //glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    //glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
    //glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
    //glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    //glEnable(GL_POLYGON_SMOOTH);
    //glEnable(GL_COLOR_MATERIAL);
    //glColorMaterial(GL_FRONT, GL_DIFFUSE);
    //glColorMaterial(GL_BACK, GL_DIFFUSE);
    //glEnable(GL_MULTISAMPLE);
    //glEnable(GL_LINE_SMOOTH);
    
    //glDisable(GL_DITHER);
    glEnable (GL_MULTISAMPLE_ARB);
    glHint (GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
    glEnable(GL_TEXTURE_2D);
	
    dl = glGenLists(4);
    
    glNewList(dl, GL_COMPILE);
    for(i=0; i < state->noOfStructs; i++){
        if(i%3 == 1){
            mat_specular[0] = 0.0;
        }
        else if(i%3 == 2){
            mat_specular[1] = 0.0;
			
        }
        else{
            mat_specular[2] = 0.0;
        }
        
        
        
        if(i%3 == 1){
            glColor3f(1.0, 0.0, 0.0);
        }
        else if(i%3 == 2){
            glColor3f(0.0, 1.0, 0.0);
            
        }
        else{
            glColor3f(0.0, 1.0, 0.0);
        }
        
        
        
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
        glMaterialfv(GL_BACK, GL_SPECULAR, mat_specular);
        DrawProt(i);
    }
    glEndList();
    
    //float r, g, b = 1.0;
    
    glNewList(dl+1, GL_COMPILE);
    for(i=0; i < state->noOfStructs; i++){
        if(i%3 == 1){
            glColor3f(1.0, 0.0, 0.0);
        }
        else if(i%3 == 2){
            glColor3f(0.0, 1.0, 0.0);
            
        }
        else{
            glColor3f(0.0, 1.0, 0.0);
        }
        DrawTrace(i);
    }
    glEndList();
    
    //float r, g, b = 1.0;
	
    glNewList(dl+2, GL_COMPILE);
    for(i=0; i < state->noOfStructs; i++){
        if(i%3 == 1){
            glColor3f(1.0, 0.0, 0.0);
        }
        else if(i%3 == 2){
            glColor3f(0.0, 1.0, 0.0);
            
        }
        else{
            glColor3f(0.0, 1.0, 0.0);
        }
        DrawRibbon(i);
    }
    glEndList();
    
    //float r, g, b = 1.0;
	
    glNewList(dl+3, GL_COMPILE);
    for(i=0; i < state->noOfStructs; i++){
        if(i%3 == 1){
            glColor3f(1.0, 0.0, 0.0);
        }
        else if(i%3 == 2){
            glColor3f(0.0, 1.0, 0.0);
            
        }
        else{
            glColor3f(0.0, 1.0, 0.0);
        }
        Draw3dRibbon(i);
    }
    glEndList();    
	
    //float r, g, b = 1.0;
    
    return;             // Initialization Went OK 
} 

/***************************************************************************
 DrawGLScene - draw everything
 ****************************************************************************/
GLvoid DrawGLScene(GLvoid)         // Here's Where We Do All The Drawing 
{ 
    // Clear The Screen And The Depth Buffer 
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
    
    
    glCallList(dl+drawmode);
    //glutSwapBuffers();
    return;           // Everything Went OK 
} 


/***************************************************************************
 ReSizeGLScene - Resize And Initialize The GL Window
 ****************************************************************************/
GLvoid ReSizeGLScene(int width, int height) 
{ 
    if (height==0)           // Prevent A Divide By Zero By 
    { 
        height=1;            // Making Height Equal One 
    } 
    glViewport(0, 0, width, height);      // Reset The Current Viewport 
    glMatrixMode(GL_PROJECTION);       // Select The Projection Matrix 
    glLoadIdentity();          // Reset The Projection Matrix 
	// Calculate The Aspect Ratio Of The Window 
    //gluPerspective(45.0f,(GLfloat)width/(GLfloat)height,0.1f,100.0f);
    //glOrtho(-1/sdepth*kWindowWidth, 1/sdepth*kWindowWidth, -1/sdepth*kWindowHeight, 1/sdepth*kWindowHeight, zNear, zFar);
    glMatrixMode(GL_MODELVIEW);       // Select The Modelview Matrix 
    glLoadIdentity();          // Reset The Modelview Matrix 
} 


/*----------------------------------------------------------------------*/
/* These functions implement a simple trackball-like motion control.    */
/*----------------------------------------------------------------------*/
/***************************************************************************
 callback_display
 ****************************************************************************/
void callback_display(void)
{
    //GLint viewport[4];
    //glGetIntegerv(GL_VIEWPORT, viewport);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	
    //gluPerspective(15.0, aspect, zNear, zFar);
	
    //glOrtho(-1/sdepth*kWindowWidth, 1/sdepth*kWindowWidth, -1/sdepth*kWindowHeight, 1/sdepth*kWindowHeight, zNear, zFar);
    glTranslatef(0.0,0.0,-100);
	
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity(); 
    //glTranslatef(0.0,0.0,-sdepth);
    glRotatef(-stheta, 1.0, 0.0, 0.0);
    glRotatef(sphi, 0.0, 0.0, 1.0);
    
	
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    DrawGLScene();
    
    //glTranslatef(0.0,0.0, sdepth);
    glutSwapBuffers();
}

/***************************************************************************
 motion
 ****************************************************************************/
void motion(int x, int y)
{
    if (leftButton)
    {
        sphi += (float)(x - downX) / 4.0;
        stheta += (float)(downY - y) / 4.0;
    }
    if (middleButton)
    {
        sdepth += (float)(downY - y) / 10.0;
    }
    downX = x;
    downY = y;
    glutPostRedisplay();
}

/***************************************************************************
 mouse
 ****************************************************************************/
void mouse(int button, int state, int x, int y)
{
    downX = x;
    downY = y;
    leftButton = ((button == GLUT_LEFT_BUTTON) && 
                  (state == GLUT_DOWN));
    middleButton = ((button == GLUT_MIDDLE_BUTTON) && 
                    (state == GLUT_DOWN));
}

/***************************************************************************
 setShader
 ****************************************************************************/
void setShader(){
    GLhandleARB my_program;
    GLhandleARB my_vertex_shader;
    GLhandleARB my_fragment_shader;
    GLsizei* elength = E_MALLOC(sizeof(GLint));
    GLchar* infolog = E_MALLOC(5000*sizeof(GLchar));
    
    // Create Shader And Program Objects
    my_program = glCreateProgramObjectARB();
    my_vertex_shader = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
    my_fragment_shader = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
    
    // Load Shader Sources
    GLcharARB* vertsrc;
    GLcharARB* fragsrc;
    switch(shademode){
        case GOOCH_SHADER:
            vertsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/gooch.vert");
            fragsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/gooch.frag");
            break;
        case PHONG_SHADER:
            vertsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/phong-use-diffuse.vert");
            fragsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/phong-use-diffuse.frag");
            break;
        case LATTICE_SHADER:
            vertsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/CH11-lattice.vert.txt");
            fragsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/CH11-lattice.frag.txt");
            break;
        case PPL2_SHADER:
            vertsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/per-pixel-lighting.vert");
            fragsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/per-pixel-lighting.frag");
            break;
        case TEST_SHADER:
            vertsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/CH14-chromaticAb.vert.txt");
            fragsrc = readShaderSrc("/Users/tom/Desktop/all_shaders/CH14-chromaticAb.frag.txt");
            break;
    }
    
    glShaderSourceARB(my_vertex_shader, 1, &vertsrc, NULL);
    glShaderSourceARB(my_fragment_shader, 1, &fragsrc, NULL);
    
    
    // Compile The Shaders
    glCompileShaderARB(my_vertex_shader);
    
    GLint blen = 0;	
    GLint slen = 0;
    GLchar* compiler_log;
    
    glGetObjectParameterivARB(my_vertex_shader, GL_OBJECT_INFO_LOG_LENGTH_ARB , &blen);
    if (blen > 1)
    {
        if ((compiler_log = (GLcharARB*)E_MALLOC(blen)) == NULL) 
        {
            printf("OOM\n");
        }
        glGetInfoLogARB(my_vertex_shader, blen, &slen, compiler_log);
        if (compiler_log!=0) {
            printf("compiler_log: %s\n", compiler_log); 
            free(compiler_log);
        }
    }
    
    
    glCompileShaderARB(my_fragment_shader);
    
    glGetObjectParameterivARB(my_fragment_shader, GL_OBJECT_INFO_LOG_LENGTH_ARB , &blen);
    if (blen > 1)
    {
        if ((compiler_log = (GLcharARB*)E_MALLOC(blen)) == NULL) 
        {
            printf("OOM\n");
        }
        glGetInfoLogARB(my_fragment_shader, blen, &slen, compiler_log);
        if (compiler_log!=0) {
            printf("compiler_log: %s\n", compiler_log); 
            free(compiler_log);
        }
    }
    
    
    // Attach The Shader Objects To The Program Object
    glAttachObjectARB(my_program, my_vertex_shader);
    glAttachObjectARB(my_program, my_fragment_shader);
    
    glGetObjectParameterivARB(my_program, GL_OBJECT_INFO_LOG_LENGTH_ARB , elength);
    if (elength > 0){
        glGetInfoLogARB(my_program, 500, elength, infolog);
        printf("%s\n", infolog);
    }
    
    
    // Link The Program Object
    glLinkProgramARB(my_program);
	
    glGetObjectParameterivARB(my_program, GL_OBJECT_INFO_LOG_LENGTH_ARB , elength);
    if (elength > 0){
        glGetInfoLogARB(my_program, 500, elength, infolog);
        printf("%s\n", infolog);
    }
    
    
    glValidateProgramARB(my_program);
    if(GL_VALIDATE_STATUS){
        // Use The Program Object Instead Of Fixed Function OpenGL
        glUseProgramObjectARB(my_program);
        
        if(shademode == GOOCH_SHADER){
			// set gooch shader parameters
            GLint loc;
            GLfloat vs[]={0.75, 0.75, 0.75};
            loc = glGetUniformLocationARB(my_program, "SurfaceColor");
            glUniform3fv(loc, 1, vs);
            GLfloat vw[]={0.6, 0.6, 0.0};
            loc = glGetUniformLocationARB(my_program, "WarmColor");
            glUniform3fv(loc, 1, vw);
            GLfloat vc[]={0.0, 0.0, 0.6};
            loc = glGetUniformLocationARB(my_program, "CoolColor");
            glUniform3fv(loc, 1, vc);
            loc = glGetUniformLocationARB(my_program, "DiffuseWarm");
            glUniform1f(loc, 0.45);
            loc = glGetUniformLocationARB(my_program, "DiffuseCool");
            glUniform1f(loc, 0.40);
        }
        if(shademode == LATTICE_SHADER){
            // set gooch shader parameters
            GLint loc;
            GLfloat vs[]={10.00, 100.0, 0.0};
            loc = glGetUniformLocationARB(my_program, "LightPosition");
            glUniform3fv(loc, 1, vs);
            GLfloat vw[]={0.6, 0.6, 0.0};
            loc = glGetUniformLocationARB(my_program, "LightColor");
            glUniform3fv(loc, 1, vw);
            GLfloat vc[]={0.0, 0.0, 100.0};
            loc = glGetUniformLocationARB(my_program, "EyePosition");
            glUniform3fv(loc, 1, vc);
            GLfloat vd[]={0.5, 0.5, 0.5};
            loc = glGetUniformLocationARB(my_program, "Specular");
            glUniform3fv(loc, 1, vd);
            GLfloat ve[]={1.0, 1.0, 1.0};
            loc = glGetUniformLocationARB(my_program, "Ambient");
            glUniform3fv(loc, 1, ve);
            loc = glGetUniformLocationARB(my_program, "Kd");
            glUniform1f(loc, 0.40);
            GLfloat vf[]={20.0, 20.0};
            loc = glGetUniformLocationARB(my_program, "Scale");
            glUniform2fv(loc, 1, vf);
            GLfloat vg[]={0.3, 0.3};
            loc = glGetUniformLocationARB(my_program, "Threshold");
            glUniform2fv(loc, 1, vg);
            GLfloat vh[]={0.5, 0.5, 0.5};
            loc = glGetUniformLocationARB(my_program, "SurfaceColor");
            glUniform3fv(loc, 1, vh);
        }
        
    }
    
    glGetObjectParameterivARB(my_program, GL_OBJECT_INFO_LOG_LENGTH_ARB , elength);
    if (elength > 0){
        glGetInfoLogARB(my_program, 500, elength, infolog);
        printf("%s\n", infolog);
    }
}

/***************************************************************************
 keyboard
 ****************************************************************************/
void keyboard(unsigned char key, int x, int y)
{
    if(key == 'm'){
        drawmode++;
        drawmode = drawmode % MAX_DRAW_MODES;
    }
    else if(key == 's'){
        shademode++;
        shademode = shademode % MAX_SHADE_MODES;
        setShader();
    }
    callback_display();
    return;
}

/***************************************************************************
 initState
 ****************************************************************************/
struct speckstate* initState(int n)
{
    struct speckstate* state = E_MALLOC(sizeof(struct speckstate));
    state->noOfStructs = n;
    state->coords = E_MALLOC(n * sizeof(struct coord*));
    return(state);
}

/* -----------------------------------------------------------------------------
 Generate a preview for file
 
 This function's job is to create preview for designated file
 ----------------------------------------------------------------------------- */ 

OSStatus GeneratePreviewForURL(void *thisInterface, QLPreviewRequestRef preview, CFURLRef url, CFStringRef contentTypeUTI, CFDictionaryRef options)
{
	state = initState(1);
    //printf("displaying %d structures \n", state->noOfStructs);
	CFStringRef filepath = CFURLGetString(url);
	state->coords[0] = pdb_read(CFStringGetCStringPtr(filepath, 0), "1xxx", '_');

    //glutInit(&argc, argv);
    //glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
    //glutInitWindowSize (kWindowWidth, kWindowHeight);
    //glutInitWindowPosition (100, 100);
    //glutCreateWindow ("WURST-o-Visison");
	
    InitGL();
	
    //setShader();
    
    glutDisplayFunc(callback_display);
    glutReshapeFunc(ReSizeGLScene);
    
    glutMouseFunc(mouse);
    glutKeyboardFunc(keyboard);
    glutMotionFunc(motion);
    
	//    glutMouseFunc(callback_mouse);
	//    glutMotionFunc(callback_motion);
    //glutIdleFunc(DrawGLScene);
    glutMainLoop();
    return noErr;
}

void CancelPreviewGeneration(void* thisInterface, QLPreviewRequestRef preview)
{
    // implement only if supported
}