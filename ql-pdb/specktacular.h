/*
 *  specktacular.h
 *  gltest
 *
 *  Created by Thomas Margraf on 21/10/2007.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#import "coord_i.h"

struct speckstate{
    char** structurenames;
    struct coord **coords;
    int* colors;
    int* highlights;
    int noOfHighlights;
    int noOfStructs;
    int* representations;
} speckstate;