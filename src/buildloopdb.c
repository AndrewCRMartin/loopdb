/************************************************************************/
/**

   \file       buildloopdb.c
   
   \version    V1.0
   \date       16.07.15
   \brief      Build a database of CDR-H3 like loops
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 2015
   \author     Dr. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College London,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is NOT IN THE PUBLIC DOMAIN, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   Reads a directory of PDB files and identifies stretches that match
   the takeoff region distances for CDR-H3 loops (i.e. H92-H94 with
   H103-H105). The mean and standard deviation distances are stored in
   distances.h which is built automatically from a directory of PDB
   files. A table containing distance ranges may be used to override
   these defaults. Output is a file containing the PDB code, residue range
   loop length (residues between the takeoff regions) and the 9 distances.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   16.07.15  Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <time.h>
#include <math.h>

#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF                160
#define MAX_CA_CA_DISTANCE_SQ   16.0  /* max CA-CA distance of 4.0      */

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *minLength, int *maxLength, BOOL *isDirectory,
                  char *distTable);
void Usage(void);
void RunAnalysis(FILE *out, PDB *pdb, int minLength, int maxLength, 
                 char *pdbCode, REAL minTable[3][3], REAL maxTable[3][3]);
void PrintResults(FILE *out, char *pdbCode, int separation, PDB *p[3], 
                  PDB *q[3], REAL distMat[3][3]);
void ProcessFile(FILE *in, FILE *out, int minLength, int maxLength,
                 char *pdbCode, REAL minTable[3][3], REAL maxTable[3][3]);
void ProcessAllFiles(FILE *out, char *dirName, int minLength, 
                     int maxLength, REAL minTable[3][3], 
                     REAL maxTable[3][3]);
void PrintHeader(FILE *out, char *dirName);
void ReadDistanceTable(char *distTable, REAL minTable[3][3],
                       REAL maxTable[3][3]);
void SetUpMinMaxTables(REAL minTable[3][3], REAL maxTable[3][3]);
BOOL ChainIsIntact(PDB *start, PDB *end);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
-  14.07.15 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        distTable[MAXBUFF];
   FILE *in         = stdin,
        *out        = stdout;
   int  minLength   = 0,
        maxLength   = 0,
        retval      = 0;
   BOOL isDirectory = FALSE;
   REAL minTable[3][3],
        maxTable[3][3];

   /* Default distance ranges for CDR-H3                                */
   SetUpMinMaxTables(minTable, maxTable);

   if(!ParseCmdLine(argc, argv, infile, outfile, &minLength, &maxLength,
                    &isDirectory, distTable))
   {
      Usage();
      return(0);
   }
   else
   {
      if(distTable[0])
      {
         ReadDistanceTable(distTable, minTable, maxTable);
      }

      if(isDirectory)
      {
         if(blOpenStdFiles(NULL, outfile, NULL, &out))
         {
            PrintHeader(out, infile);
            ProcessAllFiles(out, infile, minLength, maxLength, 
                            minTable, maxTable);
         }
      }
      else
      {
         if(blOpenStdFiles(infile, outfile, &in, &out))
         {
            char *pdbCode;
            pdbCode = blFNam2PDB(infile);
            ProcessFile(in, out, minLength, maxLength, pdbCode, 
                        minTable, maxTable);
         }
         else
         {
            return(1);
         }
      }
      
   }
         
   return(retval);
}

/************************************************************************/
/*>void PrintHeader(FILE *out, char *dirName)
   ------------------------------------------
   \param[in]   *out      Output file pointer
   \param[in]   *dirName  Name of directory being processed

   Prints a short header for the database file
*//**
-  14.07.15 Original   By: ACRM
*/
void PrintHeader(FILE *out, char *dirName)
{
   time_t tm;
   
   time(&tm);
   fprintf(out,"#PDBDIR: %s\n",dirName);
   fprintf(out,"#DATE:   %s\n",ctime(&tm));
}

   
/************************************************************************/
/*>void ProcessAllFiles(FILE *out, char *dirName, int minLength, 
                        int maxLength, REAL minTable[3][3], 
                        REAL maxTable[3][3])
   --------------------------------------------------------------
*//**
   \param[in]   *out       Output file pointer
   \param[in]   *dirName   Directory being processed
   \param[in]   minLength  Minimum loop length (0 = no limit)
   \param[in]   maxLength  Maximum loop length (0 = no limit)
   \param[in]   minTable   table of minimum distances
   \param[in]   maxTable   table of maximum distances

   Steps through all files in the specified directory and processes them
   via calls to ProcessFile()

-  14.07.15 Original   By: ACRM
*/
void ProcessAllFiles(FILE *out, char *dirName, int minLength, 
                     int maxLength, REAL minTable[3][3], 
                     REAL maxTable[3][3])
{
   DIR           *dp;
   struct dirent *dent;
   char          filename[MAXBUFF];
   FILE          *in;
   
   if((dp=opendir(dirName))!=NULL)
   {
      while((dent=readdir(dp))!=NULL)
      {
         if(dent->d_name[0] != '.')
         {
            sprintf(filename,"%s/%s",dirName,dent->d_name);
            if((in=fopen(filename, "r"))!=NULL)
            {
               char *pdbCode;
               pdbCode = blFNam2PDB(filename);
               fprintf(stderr,"%s\n", filename);
               ProcessFile(in, out, minLength, maxLength, pdbCode, 
                           minTable, maxTable);
               fclose(in);
            }
         }
      }
      closedir(dp);
   }
}

/************************************************************************/
/*>void ProcessFile(FILE *in, FILE *out, int minLength, int maxLength, 
                    char *pdbCode, REAL minTable[3][3], 
                    REAL maxTable[3][3])
   --------------------------------------------------------------------
*//**
   \param[in]   *in        Input file pointer (for PDB file)
   \param[in]   *out       Output file pointer
   \param[in]   minLength  Minimum loop length
   \param[in]   maxLength  Maximum loop length
   \param[in]   pdbCode    PDB code for this file
   \param[in]   minTable   table of minimum distances
   \param[in]   maxTable   table of maximum distances

   Obtains the PDB data and calls RunAnalysis() to do the real work

-  14.07.15 Original   By: ACRM
*/
void ProcessFile(FILE *in, FILE *out, int minLength, int maxLength, 
                 char *pdbCode, REAL minTable[3][3], REAL maxTable[3][3])
{
   PDB *pdb;
   int natoms;

   if((pdb = blReadPDBAtoms(in, &natoms))!=NULL)
   {
      /* Extract the CAs                                                */
      if((pdb = blSelectCaPDB(pdb))!=NULL)
      {
         /* Run the analysis                                            */
         RunAnalysis(out, pdb, minLength, maxLength, pdbCode, 
                     minTable, maxTable);
      
      }
      FREELIST(pdb, PDB);
   }
}



/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     int *minLength, int *maxLength, BOOL *isDirectory,
                     char *distTable)
   ---------------------------------------------------------------------
*//**
   \param[in]   argc              Argument count
   \param[in]   **argv            Argument array
   \param[out]  *infile           Input filename (or blank string)  
   \param[out]  *outfile          Output filename (or blank string) 
   \param[out]  *minLength        miniumum loop length to display   
   \param[out]  *maxLength        maxiumum loop length to display   
   \param[out]  *isDirectory      Input 'filename' is a directory   
   \param[out]  *distTable        Distance table filename           
   \return                        Success

   Parse the command line

-  14.07.15 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *minLength, int *maxLength, BOOL *isDirectory,
                  char *distTable)
{
   BOOL gotArg = FALSE;
   
   argc--;
   argv++;
   
   infile[0]    = outfile[0] = '\0';
   *minLength   = *maxLength = 0;
   *isDirectory = TRUE;
   distTable[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'm':
            argv++;
            argc--;
            if(!argc || !sscanf(argv[0], "%d", minLength))
               return(FALSE);
            break;
         case 'x':
            argv++;
            argc--;
            if(!argc || !sscanf(argv[0], "%d", maxLength))
               return(FALSE);
            break;
         case 't':
            argv++;
            argc--;
            if(!argc)
               return(FALSE);
            strncpy(distTable, argv[0], MAXBUFF);
            break;
         case 'p':
            *isDirectory = FALSE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);

         gotArg = TRUE;
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }

      argc--;
      argv++;
   }

   /* If it's a directory then we MUST have a directory name            */
   if(*isDirectory && !gotArg)
      return(FALSE);
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**
   Prints a usage message 

-  14.07.15 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nbuildloopdb V1.0 (c) 2015 UCL, Dr. Andrew C.R. \
Martin.\n");

   fprintf(stderr,"\nUsage: buildloopdb [-m minLength][-x maxLength]\
[-t disttable]\n");
   fprintf(stderr,"                   pdbdir [out.db]\n");
   fprintf(stderr,"--or--\n");
   fprintf(stderr,"       buildloopdb -p [-m minLength][-x maxLength]\
[-t disttable]\n");
   fprintf(stderr,"                   [in.pdb [out.db]]\n");
   

   fprintf(stderr,"\n                   -p Argument is a PDB file\n");
   fprintf(stderr,"                   -m Set minimum loop length\n");
   fprintf(stderr,"                   -x Set maximum loop length\n");
   fprintf(stderr,"                   -t Specify a distance table\n");

   fprintf(stderr,"\nReads a directory of PDB files and identifies \
stretches that match\n");
   fprintf(stderr,"the takeoff region distances for CDR-H3 loops (i.e. \
H92-H94 with\n");
   fprintf(stderr,"H103-H105). The mean and standard deviation distances \
are stored in\n");
   fprintf(stderr,"distances.h which is built automatically from a \
directory of PDB\n");
   fprintf(stderr,"files. Output is a file containing the PDB code, \
residue range\n");
   fprintf(stderr,"loop length (residues between the takeoff regions) \
and the 9 distances.\n");
   fprintf(stderr,"-t allows the default distance ranges to be \
overridden; the distance file\n");
   fprintf(stderr,"contains nine min/max distance pairs representing \
n0-c0, n0-c1, n0-c2,\n");
   fprintf(stderr,"n1-c0, n1-c1, n1-c2, n2-c0, n2-c1, n2-c2\n");

   fprintf(stderr,"\n-p is primarilly for testing - it builds a database \
from a single PDB\n\n");
   fprintf(stderr,"file instead of a directory of PDB files\n");

   fprintf(stderr,"\nInput/output is to standard input/output if files \
are not specified.\n");
   fprintf(stderr,"However without the -p flag, a directory name is not \
optional.\n\n");
}

/************************************************************************/
/*>void RunAnalysis(FILE *out, PDB *pdb, int minLength, int maxLength, 
                    char *pdbCode, REAL minTable[3][3], 
                    REAL maxTable[3][3])
   --------------------------------------------------------------------
*//**
   \param[in]   *out        Output file pointer
   \param[in]   *pdb        Pointer to PDB linked list
   \param[in]   minLength   Minimum loop length       
   \param[in]   maxLength   Maximum loop length       
   \param[in]   *pdbCode    PDB code for this file    
   \param[in]   minTable    table of minimum distances
   \param[in]   maxTable    table of maximum distances

   Does the real work of analyzing a structure. Steps through N-ter and
   C-ter triplets of residues to find those that match the requirements
   of the minTable and maxTable distance matrices as well as any spcified
   loop length requirements.

-  14.07.15 Original   By: ACRM
*/
void RunAnalysis(FILE *out, PDB *pdb, int minLength, int maxLength, 
                 char *pdbCode, REAL minTable[3][3], REAL maxTable[3][3])
{
   PDB  *n[3], *c[3],
        *chain,
        *nextChain;
   REAL distMat[3][3];
   int  i, j,
        separation;
   
   for(chain=pdb; chain!=NULL; chain=nextChain)
   {
      nextChain = blFindNextChainPDB(chain);
      
      /* Find an N-terminal residue                                     */
      for(n[0]=chain; n[0]!=NULL; NEXT(n[0]))
      {
         /* And find the next two                                       */
         n[1] = (n[0])?n[0]->next:NULL;
         n[2] = (n[1])?n[1]->next:NULL;

         /* If all three are valid                                      */
         if((n[2] != NULL) && (n[2]->next != NULL))
         {
            separation = 0;
            
            /* Find a C-terminal residue                                */
            for(c[0]=n[2]->next->next; c[0]!=NULL; NEXT(c[0]))
            {
               /* If the spacing between N and Cter is too long or not
                  long enough, break out
               */
               separation++;
               if(maxLength && (separation > maxLength))
                  break;

               if(separation >= minLength)
               {
                  /* And find the next two                              */
                  c[1] = (c[0])?c[0]->next:NULL;
                  c[2] = (c[1])?c[1]->next:NULL;

                  /* If all three are valid                             */
                  if((c[1] != NULL) && (c[2] != NULL))
                  {
                     BOOL badDistance = FALSE;
                  
                     if(ChainIsIntact((PDB*)(n[0]), c[2]->next))
                     {
                        /* Create the distance matrix                   */
                        for(i=0; i<3; i++)
                        {
                           for(j=0; j<3; j++)
                           {
                              distMat[i][j] = DIST(n[i], c[j]);
                              
                              if((distMat[i][j] < minTable[i][j]) ||
                                 (distMat[i][j] > maxTable[i][j]))
                              {
                                 badDistance = TRUE;
                                 break;
                              }
                           }
                        }

                        if(!badDistance)
                           PrintResults(out, pdbCode, separation, n, c, 
                                        distMat);
                     }
                  }
               }
            }  
         }
      }
   }
}


/************************************************************************/
/*>BOOL ChainIsIntact(PDB *start, PDB *end)
   ----------------------------------------
*//**
   \param[in]    *start  Start of region
   \param[in]    *end    End of region
   \return               Is intact?

   Checks whether a chain is intact (i.e. doesn't have any chain breaks)

-  16.07.15 Original   By: ACRM
*/
BOOL ChainIsIntact(PDB *start, PDB *end)
{
   PDB *p;
   
   for(p=start; p!=end; NEXT(p))
   {
      if((p!=NULL) && (p->next != NULL))
      {
         if(DISTSQ(p, p->next) > MAX_CA_CA_DISTANCE_SQ)
         {
            return(FALSE);
         }
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void PrintResults(FILE *out, char *pdbCode, int separation, 
                     PDB *n[3], PDB *c[3], REAL distMat[3][3]) 
   ------------------------------------------------------------
*//**
   \param[in]   *out          Output file pointer
   \param[in]   *pdbCode      PDB code
   \param[in]   separation    loop length
   \param[in]   *n[]          N-ter three PDB pointers
   \param[in]   *c[]          C-ter three PDB pointers
   \param[in]   *distMat[][]  Distance matrix

   Prints the results for a loop already determined to match criteria

-  14.07.15 Original   By: ACRM
*/
void PrintResults(FILE *out, char *pdbCode, int separation, 
                  PDB *n[3], PDB *c[3], REAL distMat[3][3]) 
{
   char resid1[16],
        resid2[16];
   int  i, j;
               
   MAKERESID(resid1, n[0]);
   MAKERESID(resid2, c[2]);
   
   fprintf(out, "%s %s %s %d ", (pdbCode!=NULL)?pdbCode:"",
           resid1, resid2, separation);
   for(i=0; i<3; i++)
   {
      for(j=0; j<3; j++)
      {
         fprintf(out, "%.3f ", distMat[i][j]);
      }
   }
   fprintf(out, "\n");
}


/************************************************************************/
/*>void ReadDistanceTable(char *distTable, REAL minTable[3][3], 
                          REAL maxTable[3][3])
   ------------------------------------------------------------
*//**
   \param[in]   *distTable    Distance table filename
   \param[out]  minTable[][]  table of minimum distances
   \param[out]  maxTable[][]  table of maximum distances

   Reads a user-specified distance matrix table instead of using the
   defaults coded in distances.h

-  14.07.15 Original   By: ACRM
*/
void ReadDistanceTable(char *distTable, REAL minTable[3][3],
                       REAL maxTable[3][3])
{
   FILE *fp = NULL;
   char buffer[MAXBUFF];
   int  i = 0,
        j = 0;
   
   if((fp=fopen(distTable, "r"))!=NULL)
   {
      while(fgets(buffer, MAXBUFF, fp))
      {
         char *chp;
         TERMINATE(buffer);
         if((chp = strchr(buffer, '#'))!=NULL)
            *chp = '\0';
         KILLTRAILSPACES(buffer);
         if(strlen(buffer))
         {
            sscanf(buffer,"%lf %lf", 
                   &(minTable[i][j]), &(maxTable[i][j]));
            if((++j)==3)
            {
               j=0;
               i++;
            }
         }
      }
      
      fclose(fp);
   }
}


/************************************************************************/
/*>void SetUpMinMaxTables(REAL minTable[3][3], REAL maxTable[3][3])
   ----------------------------------------------------------------
*//**
   \param[out]  minTable[][]  table of minimum distances
   \param[out]  maxTable[][]  table of maximum distances

   Initializes the minimum and maximum distance matrices based on
   mean and standard deviation distances in distances.h

-  15.07.15 Original   By: ACRM
*/
void SetUpMinMaxTables(REAL minTable[3][3], REAL maxTable[3][3])
{
   int i, j;

   #include "distances.h"

   for(i=0; i<3; i++)
   {
      for(j=0; j<3; j++)
      {
         minTable[i][j] = means[i][j] - sdMult * sds[i][j];
         maxTable[i][j] = means[i][j] + sdMult * sds[i][j];
      }
   }
}
