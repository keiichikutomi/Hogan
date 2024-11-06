/*MODEL001.C FOR WIN32 SINCE 2015.02.28.FUKUSHIMA.*/
/*LAST CHANGE:07-Mar-2015.*/
/*DXF-TO-INP CONVERSION*/
/*SUPPORTED ENTITIES: LINE, 3DFACE, HATCH*/

// TODO: SUPPORT POLYLINE


enum ENTITYTYPE {
  UNKNOWN,
  LINE,
  TDFACE,
  HATCH,
};

char *extension(const char *filename);
int inputorganizationdxf(FILE *fin,
                         struct organ *org,
                         struct viewparam *vp);
struct onode *coordnode(struct organ *org,double x,double y,double z,double eps);
struct osect *searchsection(struct organ *org,int code);
struct onode *nodeinbox(struct organ *org,struct onode n1,struct onode n2,double eps,int *size);
double *direction(struct onode n1,struct onode n2);
double distance(struct onode n1,struct onode n2);
char parallel(double *vec1,double *vec2,double eps);
void sortnode(struct onode *before,double *keys,int size);
struct onode *onnode(struct organ *org,struct oelem elem,double eps,int *size);
// DIVIDE
struct oelem *divideatnode(struct organ *org,struct oelem *elem,struct onode node,int position);
struct oelem *divideatonnode(struct organ *org,struct oelem elem,double eps,int *size);
struct oelem *cutbyelem(struct organ *org,struct oelem cutter,struct oelem cuttee,char cross,double eps);
void intersectall(struct organ *org,double eps);


char *extension(const char *filename)
/*RETURN EXTENSION OF FILENAME WITHOUT ".".*/
{
  char *rtn;
  char *tmp = strtok(strdup(filename), ".");
  while (tmp != NULL)
  {
    rtn = tmp;
    tmp = strtok(NULL, ".");
  }
  return rtn;
} /*extension*/

int inputorganizationdxf(FILE *fin,
                         struct organ *org,
                         struct viewparam *vp,
                         double factor,
                         double eps)
/*INPUT ORGAN FROM DXFFILE.*/
/*FACTOR: UNIT LENGTH OF INP FILE DIVIDED BY THAT OF DXF FILE.(e.g. 0.001: [mm]->[m])*/
/*EPS   : COORDINATE TOLERANCE.[m]*/
//DXF GROUPE CODE
//  COMMON
//    0:ENTITY TYPE
//    8:LAYERNAME(ETYPE+SECTCODE)
//  LINE
//    10:X1, 20:Y1, 30:Z1
//    11:X2, 21:Y2, 31:Z2
//  3DFACE
//    10:X1, 20:Y1, 30:Z1
//    11:X2, 21:Y2, 31:Z2
//    12:X3, 22:Y3, 32:Z3
//    13:X4, 23:Y4, 33:Z4
//  HATCH
//    10:- , 20:- , 30:Z
//    93:NUMBER OF VERTICES(=ENODS)
//    10:X1, 20:Y1
//    10:X2, 20:Y2
//    10:X3, 20:Y3
//    10:X4, 20:Y4
{
  char str[256];
  char data[80];
  char num[3];
  int i,j;
  char codeline,endsec,add;
  int groupcode,scode;
  int ind;
  int nelem,nline;
  double *cx,*cy,*cz;
  long entstart;
  int size;
  int hedge,hcoord;
  struct onode *n;
  struct oelem *elems,*pelem;
  struct osect *ps;
  ENTITYTYPE enttype=UNKNOWN;

  org->code=0;
  org->loff=0;

  org->nnode=0; /*INITIALIZATION.*/
  org->nelem=0;
  org->nprop=0;
  org->nsect=0;
  //org->npile=0;

  org->ntext=0;
  org->texts=NULL;

  org->ai.Cox=0.2;
  org->ai.Coy=0.2;
  org->ai.Z =1.0;
  org->ai.T1=0.02;
  org->ai.Tc=0.6;
  org->ai.nfloor=19;
  org->ai.lbound=(double *)malloc((org->ai.nfloor+1)*sizeof(double));
  org->ai.nfloor=-1;


  fseek(fin,0L,SEEK_SET);
  while(1) // SEARCH FOR ENTITIES SECTION
  {
    if(fgets(str,256,fin)==NULL)
    {
      return 1;
    }
    if(!strncmp(str,"SECTION",7))
    {
      fgets(str,256,fin);
      fgets(str,256,fin);
      if(!strncmp(str,"ENTITIES",8))
      {
        break;
      }
    }
  }
  entstart=ftell(fin);
  codeline=1;
  groupcode=0;
  endsec=0;
  size=0;
  while(fgets(str,256,fin)!=NULL && !endsec) // COUNT MAXIMUM NUMBER OF ELEMENTS
  {
    if(codeline)
    {
      groupcode=(int)strtol(str,NULL,10);
      codeline--;
    }
    else
    {
      switch(groupcode)
      {
        case 0:
          if(!strncmp(str,"ENDSEC",6))
          {
            endsec=1;
          }
          else
          {
            size++;
          }
          break;
      }
      codeline++;
    }
  }
  fseek(fin,entstart,SEEK_SET);
  codeline=1;
  groupcode=0;
  endsec=0;
  add=0;
  nelem=0;
  nline=0;
  cx=(double *)malloc(4*sizeof(double));
  cy=(double *)malloc(4*sizeof(double));
  cz=(double *)malloc(4*sizeof(double));
  elems=(struct oelem *)malloc(size*sizeof(struct oelem));
  pelem=(struct oelem *)malloc(sizeof(struct oelem));
  while(fgets(str,256,fin)!=NULL && !endsec) // CREATE ELEMENTS
  {
    nline++;
    if(codeline)
    {
      groupcode=(int)strtol(str,NULL,10);
      codeline--;
    }
    else
    {
      switch(groupcode)
      {
        case 0:
          if((enttype==TDFACE || enttype==HATCH) && pelem->nnod<3)
          {
            add=0;
          }
          if(add)
          {
            pelem->code=1001+nelem;
            pelem->loff=nelem;
            pelem->nods=(struct onode **)malloc(pelem->nnod*sizeof(struct onode *));
            pelem->bans=(struct obans *)malloc(pelem->nban*sizeof(struct obans));
            for(i=0;i<pelem->nban;i++)
            {
              pelem->bans->code=i+1;
              pelem->bans->nnod=pelem->nnod;
              pelem->bans->nods=(struct onode **)malloc(pelem->bans->nnod*sizeof(struct onode *));
            }
            ind=0;
            for(i=0;i<pelem->nnod;i++)
            {
              if(enttype==HATCH)
              {
                n=coordnode(org,cx[i]*factor,cy[i]*factor,cz[0]*factor,eps);
              }
              else
              {
                n=coordnode(org,cx[i]*factor,cy[i]*factor,cz[i]*factor,eps);
              }
              if(i>=2 && (*(pelem->nods+i-1))->code==n->code) // ELIMINATE DUPLICATED ENOD
              {
                pelem->nnod--;
                for(j=0;j<pelem->nban;j++)
                {
                  (pelem->bans+j)->nnod--;
                }
                ind++;
              }
              else
              {
                *(pelem->nods+i-ind)=(struct onode *)
                                      malloc(sizeof(struct onode));
                **(pelem->nods+i-ind)=*n;
                for(j=0;j<pelem->nban;j++)
                {
                  *((pelem->bans+j)->nods+i-ind)=(struct onode *)
                                                  malloc(sizeof(struct onode));
                  **((pelem->bans+j)->nods+i-ind)=*n;
                }
              }
            }
            pelem->bonds=(signed char *)malloc(6*pelem->nnod*sizeof(char));
            for(i=0;i<6*(pelem->nnod);i++)
            {
              *(pelem->bonds+i)=0;
            }
            pelem->cangle=0.0;
            for(i=0;i<2;i++)
            {
              for(j=0;j<6;j++)
              {
                pelem->initial[i][j]=0.0;
              }
            }
            ps=searchsection(org,scode);
            pelem->sect=ps;
            if(pelem->nnod>=2)
            {
              *(elems+nelem)=*pelem;
              nelem++;
            }
          }

          add=0;
          enttype=UNKNOWN;
          for(i=0;i<4;i++)
          {
            cx[i]=0.0;
            cy[i]=0.0;
            cz[i]=0.0;
          }

          if(!strncmp(str,"ENDSEC",6))
          {
            endsec=1;
          }
          else if(!strncmp(str,"LINE",4))
          {
            enttype=LINE;

            pelem->code=0;
            pelem->loff=0;
            pelem->type=0;
            pelem->nnod=2;
            pelem->nban=0;
            pelem->color.line=setrgbcolor(255,255,255);

            pelem->lface[0]=0.0;
            pelem->lface[1]=0.0;
            pelem->hface[0]=0.0;
            pelem->hface[1]=0.0;
            pelem->wrect[0]=0.0;
            pelem->wrect[1]=0.0;
            add=1;
          }
          else if(!strncmp(str,"3DFACE",6))
          {
            enttype=TDFACE;

            pelem->code=0;
            pelem->loff=0;
            pelem->type=0;
            pelem->nnod=2;
            pelem->nban=1;
            pelem->color.line=setrgbcolor(255,255,255);

            pelem->lface[0]=0.0;
            pelem->lface[1]=0.0;
            pelem->hface[0]=0.0;
            pelem->hface[1]=0.0;
            pelem->wrect[0]=0.0;
            pelem->wrect[1]=0.0;
            add=1;
          }
          else if(!strncmp(str,"HATCH",5))
          {
            enttype=HATCH;
            hedge=0;
            hcoord=-1;

            pelem->code=0;
            pelem->loff=0;
            pelem->type=0;
            pelem->nnod=2;
            pelem->nban=1;
            pelem->color.line=setrgbcolor(255,255,255);

            pelem->lface[0]=0.0;
            pelem->lface[1]=0.0;
            pelem->hface[0]=0.0;
            pelem->hface[1]=0.0;
            pelem->wrect[0]=0.0;
            pelem->wrect[1]=0.0;
            add=1;
          }
          break;
        case 8:
          if(!strncmp(str,"COLUMN",6))
          {
            for(i=0;i<3;i++) num[i]=str[6+i];
            scode=(int)strtol(num,NULL,10);
            pelem->type=COLUMN;
          }
          else if(!strncmp(str,"GIRDER",6))
          {
            for(i=0;i<3;i++) num[i]=str[6+i];
            scode=(int)strtol(num,NULL,10);
            pelem->type=GIRDER;
          }
          else if(!strncmp(str,"BRACE",5))
          {
            for(i=0;i<3;i++) num[i]=str[5+i];
            scode=(int)strtol(num,NULL,10);
            pelem->type=BRACE;
          }
          else if(!strncmp(str,"SLAB",4))
          {
            for(i=0;i<3;i++) num[i]=str[4+i];
            scode=(int)strtol(num,NULL,10);
            pelem->type=SLAB;
          }
          else if(!strncmp(str,"WALL",4))
          {
            for(i=0;i<3;i++) num[i]=str[4+i];
            scode=(int)strtol(num,NULL,10);
            pelem->type=WALL;
          }
          else
          {
            add=0;
          }
          break;
        case 10:
          if(enttype==HATCH)
          {
            if(hcoord>=0 && hcoord<=3)
            {
              cx[hcoord]=strtod(str,NULL);
            }
            break;
          }
          cx[0]=strtod(str,NULL);
          break;
        case 20:
          if(enttype==HATCH)
          {
            if(hcoord>=0 && hcoord<=3)
            {
              cy[hcoord]=strtod(str,NULL);
              hcoord++;
              pelem->nnod=hcoord;
            }
            break;
          }
          cy[0]=strtod(str,NULL);
          break;
        case 30:
          cz[0]=strtod(str,NULL);
          break;
        case 11:
          cx[1]=strtod(str,NULL);
          break;
        case 21:
          cy[1]=strtod(str,NULL);
          break;
        case 31:
          cz[1]=strtod(str,NULL);
          break;
        case 12:
          cx[2]=strtod(str,NULL);
          break;
        case 22:
          cy[2]=strtod(str,NULL);
          break;
        case 32:
          cz[2]=strtod(str,NULL);
          pelem->nnod=3;
          break;
        case 13:
          cx[3]=strtod(str,NULL);
          break;
        case 23:
          cy[3]=strtod(str,NULL);
          break;
        case 33:
          cz[3]=strtod(str,NULL);
          pelem->nnod=4;
          break;
        case 93:
          hedge=strtol(str,NULL,10);
          hcoord=0;
      }
      codeline++;
    }
    currentpivot(nelem,size);
  }
  org->nelem=nelem;
  if(size!=nelem)
  {
    org->elems=(struct oelem *)malloc(nelem*sizeof(struct oelem));
    for(i=0;i<nelem;i++)
    {
      *(org->elems+i)=*(elems+i);
    }
    free(elems);
  }
  else
  {
    org->elems=elems;
  }
  // TODO: without setting elem->sect again, ''section list on/off'' seems not to work well.
  for(i=0;i<org->nelem;i++)
  {
    for(j=0;j<org->nsect;j++)
    {
      if((org->elems+i)->sect->code==(org->sects+j)->code)
      {
        (org->elems+i)->sect=(org->sects+j);
        break;
      }
    }
  }
  intersectall(org,eps);
  free(pelem);
  free(cx);
  free(cy);
  free(cz);
  return 0;
}/*inputorganizationdxf*/

struct onode *coordnode(struct organ *org,double x,double y,double z,double eps)
/*RETURN ONODE AT (x,y,z).*/
/*IF NOT EXIST, CREATE AND ADD TO ORGAN.*/
{
  long int nnode;
  int j;
  struct onode *cnode;
  struct oconf *oc;

  nnode=org->nnode;


  cnode=(struct onode *)malloc(sizeof(struct onode));

  cnode->d[0]=x;
  cnode->d[1]=y;
  cnode->d[2]=z;
//  cnode->pile=NULL;

  if(nnode==0)
  {
    oc=(struct oconf *)malloc(6*sizeof(struct oconf));
    for(j=0;j<6;j++) /*COPY CONFS*/
    {
      (oc+j)->iconf=0;
      (oc+j)->value=0.0;
    }
    org->nodes=(struct onode *)malloc(sizeof(struct onode));
    org->confs=(struct oconf *)malloc(6*sizeof(struct oconf));
    org->nnode=1;
    for(j=0;j<6;j++)
    {
      *(org->confs+j)=oc[j];
    }
    cnode->code=101;
    cnode->loff=0;
    *(org->nodes+0)=*cnode;
    return cnode;
  }

  for(j=0;j<nnode;j++) /*CANCEL IF OVERLAPPED*/
  {
    if(cnode->d[0]>(org->nodes+j)->d[0]-eps &&
       cnode->d[0]<(org->nodes+j)->d[0]+eps &&
       cnode->d[1]>(org->nodes+j)->d[1]-eps &&
       cnode->d[1]<(org->nodes+j)->d[1]+eps &&
       cnode->d[2]>(org->nodes+j)->d[2]-eps &&
       cnode->d[2]<(org->nodes+j)->d[2]+eps)
    {
      return org->nodes+j;
    }
  }
  cnode->code=(org->nodes+nnode-1)->code+1;
  cnode->loff=nnode;
  return addnode(*cnode,NULL,org);
}/*coordnode*/

struct osect *searchsection(struct organ *org,int code)
/*RETURN OSECT WITH DESIGNATED CODE.*/
/*IF NOT EXIST, CREATE AND ADD TO ORGAN.*/
{
  int i;
  long int nsect;
  struct osect *os,*ps;

  for(i=0;i<org->nsect;i++)
  {
    if((org->sects+i)->code==code)
    {
      return org->sects+i;
    }
  }

  os=(struct osect *)malloc(sizeof(struct osect));
  os->code=code;

  os->name=(char *)malloc(1*sizeof(char));
  strcpy(os->name,"\0");
  os->role=ROLENULL;
  os->type=TYPENULL;
  os->nfig=0;
  os->figs=(struct ofigs *)malloc(os->nfig*sizeof(struct ofigs));
  os->dflag=1;
  os->dcolor.r=255;
  os->dcolor.g=255;
  os->dcolor.b=255;

  os->E   =0.0;
  os->poi =0.0;
  os->area=0.0;
  os->Ixx =0.0;
  os->Iyy =0.0;
  os->Jzz =0.0;
  os->hiju[0]=0.0;
  os->hiju[1]=0.0;
  os->hiju[2]=0.0;

  os->ppc.npcurve=0;

  os->lload[0]=0.0;
  os->lload[1]=0.0;
  os->lload[2]=0.0;
  os=(struct osect *)malloc(sizeof(struct osect));
  os->code=code;

  os->name=(char *)malloc(1*sizeof(char));
  strcpy(os->name,"\0");
  os->role=ROLENULL;
  os->type=TYPENULL;
  os->nfig=0;
  os->figs=(struct ofigs *)malloc(os->nfig*sizeof(struct ofigs));
  os->dflag=1;
  os->dcolor.r=255;
  os->dcolor.g=255;
  os->dcolor.b=255;

  os->E   =0.0;
  os->poi =0.0;
  os->area=0.0;
  os->Ixx =0.0;
  os->Iyy =0.0;
  os->Jzz =0.0;
  os->hiju[0]=0.0;
  os->hiju[1]=0.0;
  os->hiju[2]=0.0;

  os->ppc.npcurve=0;

  os->lload[0]=0.0;
  os->lload[1]=0.0;
  os->lload[2]=0.0;
  ps=addsection(os,org);
  return ps;
}/*searchsection*/

struct onode *nodeinbox(struct organ *org,struct onode n1,struct onode n2,double eps,int *size)
/*RETURN ONODES IN BOX WHICH HAS n1 & n2 ON ITS CORNER*/
{
  int i,nnode;
  double *maxcoord, *mincoord;
  struct onode *nodes,*rtn;

  maxcoord=(double *)malloc(3*sizeof(double));
  mincoord=(double *)malloc(3*sizeof(double));
  for(i=0;i<3;i++)
  {
    if(n1.d[i] < n2.d[i])
    {
      mincoord[i]=n1.d[i];
      maxcoord[i]=n2.d[i];
    }
    else
    {
      mincoord[i]=n2.d[i];
      maxcoord[i]=n1.d[i];
    }
  }
  nnode=0;
  nodes=(struct onode *)malloc(org->nnode*sizeof(struct onode));
  for(i=0;i<org->nnode;i++)
  {
    if(maxcoord[0]>=(org->nodes+i)->d[0]-eps &&
       mincoord[0]<=(org->nodes+i)->d[0]+eps &&
       maxcoord[1]>=(org->nodes+i)->d[1]-eps &&
       mincoord[1]<=(org->nodes+i)->d[1]+eps &&
       maxcoord[2]>=(org->nodes+i)->d[2]-eps &&
       mincoord[2]<=(org->nodes+i)->d[2]+eps)
    {
      *(nodes+nnode)=*(org->nodes+i);
      nnode++;
    }
  }
  rtn=(struct onode *)malloc(nnode*sizeof(struct onode));
  for(i=0;i<nnode;i++)
  {
    *(rtn+i)=*(nodes+i);
  }
  free(nodes);
  free(maxcoord);
  free(mincoord);
  if(size!=NULL)
  {
    *size=nnode;
  }
  return rtn;
}/*nodeinbox*/

double *direction(struct onode n1,struct onode n2)
/*VECTOR FROM n1 TO n2*/
{
  int i;
  double *vec;
  vec=(double *)malloc(3*sizeof(double));
  for(i=0;i<3;i++)
  {
    vec[i]=n2.d[i]-n1.d[i];
  }
  return vec;
}/*direction*/

double distance(struct onode n1,struct onode n2)
/*DISTANCE BETWEEN n1 AND n2*/
{
  int i;
  double sum=0.0;
  for(i=0;i<3;i++)
  {
    sum+=(n2.d[i]-n1.d[i])*(n2.d[i]-n1.d[i]);
  }
  return sqrt(sum);
}/*distance*/

char parallel(double *vec1,double *vec2,double eps)
/*CHECK IF VEC1 & VEC2 IS PARALLEL OR NOT*/
/*1: PARALLEL*/
/*0: NOT PARALLEL*/
{
  int i,j;
  double val;
  for(i=0;i<3;i++)
  {
    j=i+1;
    if(j>=3)
    {
      j-=3;
    }
    val=vec1[i]*vec2[j]-vec1[j]*vec2[i];
    if(val>eps || val<-eps)
    {
      return 0;
    }
  }
  return 1;
}/*parallel*/

void sortnode(struct onode *before,double *keys,int size,char reverse)
/*SORT ONODE WITH KEYS BY BUBBLE SORT*/
{
  int i,j;
  double tmpk;
  struct onode tmp;

  for(i=0;i<size;i++)
  {
    for(j=size-1;j>i;j--)
    {
      if((reverse && keys[j-1]<keys[j]) ||
          (!reverse && keys[j-1]>keys[j]))
      {
        tmpk=*(keys+j);
        *(keys+j)=*(keys+j-1);
        *(keys+j-1)=tmpk;
        tmp=*(before+j);
        *(before+j)=*(before+j-1);
        *(before+j-1)=tmp;
      }
    }
  }
}

struct onode *onnode(struct organ *org,struct oelem elem,double eps,int *size)
/*NODES ON THE ELEMENT*/
{
  int i,j;
  struct onode *cand,*rtn;
  double *d,*d0,*keys,*newkeys;
  double l;
  char fr=0;
  int num=0;
  if(size==NULL)
  {
    size=(int *)malloc(sizeof(int));
    fr=1;
  }
  d0=direction(**(elem.nods+0),**(elem.nods+1));
  cand=nodeinbox(org,**(elem.nods+0),**(elem.nods+1),eps,size);
  keys=(double *)malloc((*size)*sizeof(double));
  for(i=0;i<*size;i++)
  {
    if((*(cand+i)).code==(**(elem.nods+0)).code ||
       (*(cand+i)).code==(**(elem.nods+1)).code)
    {
      *(keys+i)=-1;
      continue;
    }
    d=direction(**(elem.nods+0),*(cand+i));
    if(parallel(d,d0,eps))
    {
      l=distance(**(elem.nods+0),*(cand+i));
      *(keys+i)=l;
      num++;
    }
    else
    {
      *(keys+i)=-1;
    }
    free(d);
  }
  if(num==0)
  {
    return NULL;
  }
  rtn=(struct onode *)malloc(num*sizeof(onode));
  newkeys=(double *)malloc(num*sizeof(double));
  i=0;
  j=0;
  while(1)
  {
    if(keys[i]>0.0)
    {
      *(rtn+j)=*(cand+i);
      *(newkeys+j)=*(keys+i);
      j++;
      if(j>=num)
      {
        break;
      }
    }
    i++;
  }
  sortnode(rtn,newkeys,num,1);
  if(fr)
  {
    free(size);
  }
  else
  {
    *size=num;
  }
  free(keys);
  free(newkeys);
  free(cand);
  free(d0);
  return rtn;
}/*onnode*/

struct oelem *divideatnode(struct organ *org,struct oelem elem,struct onode node,int position)
/* DIVIDE ELEM AT NODE*/
{
  int i,j;
  struct oelem *rtn,*pelem;

  if(elem.nnod!=2)
  {
    return NULL;
  }
  for(i=0;i<elem.nnod;i++)
  {
    if((**(elem.nods+i)).code==node.code)
    {
      return NULL;
    }
  }
  rtn=(struct oelem *)malloc(2*sizeof(oelem));
  *(rtn+0)=elem;
  pelem=(struct oelem *)malloc(sizeof(struct oelem));
  pelem->code=(org->elems+org->nelem-1)->code+1;
  pelem->loff=org->nelem;
  pelem->type=elem.type;
  pelem->nnod=2;
  pelem->nban=0;
  pelem->bans=(struct obans *)malloc(pelem->nban*sizeof(struct obans));
  pelem->color.line=setrgbcolor(255,255,255);
  pelem->lface[0]=0.0;
  pelem->lface[1]=0.0;
  pelem->hface[0]=0.0;
  pelem->hface[1]=0.0;
  pelem->wrect[0]=0.0;
  pelem->wrect[1]=0.0;
  pelem->bonds=(signed char *)malloc(6*pelem->nnod*sizeof(char));
  for(i=0;i<6*(pelem->nnod);i++)
  {
    *(pelem->bonds+i)=*(elem.bonds+i);
  }
  pelem->cangle=0.0;
  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++)
    {
      pelem->initial[i][j]=0.0;
    }
  }
  pelem->sect=elem.sect;
  pelem->nods=(struct onode **)malloc(pelem->nnod*sizeof(struct onode *));
  for(i=0;i<pelem->nnod;i++)
  {
    *(pelem->nods+i)=(struct onode *)malloc(sizeof(struct onode));
  }
  switch(position)
  {
    case 1:
    case -1:
      **(pelem->nods+0)=node;
      **(pelem->nods+1)=**(elem.nods+1);
      **(elem.nods+1)=node;
      break;
    case 0:
      **(pelem->nods+0)=node;
      **(pelem->nods+1)=**(elem.nods+0);
      break;
    case 2:
      **(pelem->nods+0)=**(elem.nods+1);
      **(pelem->nods+1)=node;
      break;
  }
  *(rtn+1)=*pelem;
  addelement(pelem,org);
  free(pelem);
  return rtn;
}/*divideatnode*/

struct oelem *divideatonnode(struct organ *org,struct oelem elem,double eps,int *size)
/* DIVIDE ELEM AT NODES ON IT*/
{
  int i;
  int *nb;
  struct onode *ns;
  struct oelem *els,*rtn;

  nb=(int *)malloc(sizeof(int));

  ns=onnode(org,elem,eps,nb);
  if(ns==NULL)
  {
    free(ns);
    free(nb);
    rtn=(struct oelem *)malloc(sizeof(struct oelem));
    *(rtn+0)=elem;
    *size=1;
    return rtn;
  }
  rtn=(struct oelem *)malloc((*nb+1)*sizeof(struct oelem));
  *(rtn+0)=elem;
  for(i=0;i<*nb;i++)
  {
    els=divideatnode(org,elem,*(ns+i),1);
    *(rtn+i+1)=*(els+1);
    free(els);
  }
  *size=*nb+1;
  free(ns);
  free(nb);
  return rtn;
}/*divideatonnode*/

struct oelem *cutbyelem(struct organ *org,struct oelem cutter,struct oelem cuttee,char cross,double eps)
/*DIVIDE ELEM(cuttee) BY ANOTHER ELEM(cutter)*/
/*IF cross==1, ONLY WHEN cutter & cuttee CROSSES*/
{
  int ok;
  double k1,k2,d;
  double *d1;
  struct onode *n;
  struct oelem *rtn;
  struct line l1,l2;
  if(cutter.nnod!=2 || cuttee.nnod!=2)
  {
    return NULL;
  }
  l1.ends[0]=**(cutter.nods+0);
  l1.ends[1]=**(cutter.nods+1);
  l2.ends[0]=**(cuttee.nods+0);
  l2.ends[1]=**(cuttee.nods+1);
  //*distance=distancedotdot(n1,n2);
  ok=distancelineline(l1,l2,&k1,&k2,NULL,NULL);
  if(!ok || d>eps)
  {
    return NULL;
  }
  if(!cross || ((-eps<=k1 && k1<=1.0+eps) && (-eps<=k2 && k2<=1.0+eps)))
  {
    d1=direction(**(cutter.nods+0),**(cutter.nods+1));
    n=coordnode(org,
                (**(cutter.nods+0)).d[0]+k1*d1[0],
                (**(cutter.nods+0)).d[1]+k1*d1[1],
                (**(cutter.nods+0)).d[2]+k1*d1[2],eps);
    if(k2<-eps)
    {
      rtn=divideatnode(org,cuttee,*n,0);
    }
    else if(-eps<=k2 && k2<=1.0+eps)
    {
      rtn=divideatnode(org,cuttee,*n,1);
    }
    else
    {
      rtn=divideatnode(org,cuttee,*n,2);
    }
    free(d1);
    return rtn;
  }
  else
  {
    return NULL;
  }
}/*cutbyelem*/

void intersectall(struct organ *org,double eps)
/*DIVIDE ELEMS AT INTERSECTION*/
{
  char err[256];
  int i,j;
  int nelem,loff;
  int size;
  int total,tmp;
  struct oelem *els=NULL;
  int *checked;
  nelem=org->nelem;
  if(nelem<=1)
  {
    return;
  }
  checked=(int *)malloc(nelem*(nelem-1)*sizeof(int));
  if(checked==NULL)
  {
    return;
  }
  total=0;
  for(i=0;i<nelem;i++)
  {
    if((*(org->elems+i)).nnod==2)
    {
      loff=(*(org->elems+i)).loff;
      els=divideatonnode(org,*(org->elems+i),eps,&size);
      for(j=0;j<size;j++)
      {
        *(checked+total)=(int)(*(els+j)).loff;
        total++;
      }
      break;
    }
  }
  if(i==nelem) // NO LINE ELEM
  {
    if(els!=NULL)
    {
      free(els);
    }
    free(checked);
    return;
  }
  for(i=loff+1;i<nelem;i++)
  {
    if((*(org->elems+i)).nnod!=2)
    {
      continue;
    }
    tmp=0;
    for(j=0;j<total;j++)
    {
      els=cutbyelem(org,*(org->elems+i),*(org->elems+*(checked+j)),1,eps);
      if(els==NULL)
      {
        continue;
      }
      *(checked+total+tmp)=(*(els+1)).loff;
      tmp++;
    }
    total+=tmp;
    els=divideatonnode(org,*(org->elems+i),eps,&size);
    if(els==NULL)
    {
      continue;
    }
    for(j=0;j<size;j++)
    {
      *(checked+total)=(*(els+j)).loff;
      total++;
    }
    currentpivot(i+1,nelem);
  }
  if(els!=NULL)
  {
    free(els);
  }
  free(checked);
}/*intersectall*/
