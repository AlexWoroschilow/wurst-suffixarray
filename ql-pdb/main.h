/*
 *  main.h
 *  gltest
 *
 *  Created by Thomas Margraf on 21/6/2007.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <GLUT/glut.h> 

#define kWindowWidth 600 
#define kWindowHeight 400 

GLvoid InitGL(GLvoid);
GLvoid DrawGLScene(GLvoid);
GLvoid ReSizeGLScene(int Width, int Height);
void spinScene(void);

         