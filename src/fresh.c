#include "fargo3d.h"
#include <J_jupiter.h>

void send2cpu() {
#ifdef GPU
  Field* g;
  g = ListOfGrids;
  printf("\nCopying Fields--------------------------------------\n");
  while(g != NULL){
    if(Dev2Host3D(g)){
      printf("Error in send2cpu() with the field %s\n.", g->name);
      exit(-1);
    }
    else {
      printf("Field %s has been copied (Dev2Host)\n", g->name);
    }
    g = g->next;
  }
  printf("\n--------------------------------------------------\n\n");
#endif
}

void send2gpu() {
#ifdef GPU
  Field* g;
  g = ListOfGrids;
  printf("\nCopying Fields--------------------------------------\n");
  while(g != NULL){
    if(Host2Dev3D(g)){
      printf("Error in send2gpu() with the field %s\n.", g->name);
      exit(-1);
    }
    else {
      printf("Field %s has been copied (Host2Dev)\n", g->name);
    }
    g = g->next;
  }
  printf("\n--------------------------------------------------\n\n");
#endif
}

void Input2D_CPU(Field2D *field, int line, const char *string){
  int status;
#ifdef GPU
  if(!(*(field->fresh_cpu))) {
    status = Dev2Host2D(field);
    if(status!=0) {
      printf("Error %d in Dev2Host2D()! Field: %s\n",status,field->name);
      printf("called from line %d of file %s\n", line, string);
      exit(EXIT_FAILURE);
    }
  }
//  else {
//    printf("Field %s is up to date on Host.\n",field->name);
//  }
  *(field->fresh_cpu) = YES;
#endif
  return;
}

void Input2D_GPU(Field2D *field, int line, const char *string){
  int status;
#ifdef GPU
  if(!(*(field->fresh_gpu))) {
    status = Host2Dev2D(field);
    if(status!=0) {
      printf("Error %d in Host2Dev2D()! Field: %s level %d Pitch2D %d\n",status,field->name,field->desc->level,field->desc->Pitch2D);
      printf("called from line %d of file %s\n", line, string);
      exit(EXIT_FAILURE);
    }
  }
//  else {
//    printf("Field %s is up to date on Device.\n",field->name);
//  }
  *(field->fresh_gpu) = YES;
#endif
  return;
}

void Input2DInt_CPU(FieldInt2D *field, int line, const char *string){
  int status;
#ifdef GPU
  if(!(*(field->fresh_cpu))) {
    status = Dev2Host2DInt(field);
    //    printf("%d\n",status);
    //    printf("Copying %s from Dev to Host\n",field->name);
    if(status!=0) {
      printf("Error %d in Dev2Host2DInt()! Field: %s\n",status,field->name);
      printf("called from line %d of file %s\n", line, string);
      exit(EXIT_FAILURE);
    }
  }
//  else {
//    printf("Field %s is up to date on Host.\n",field->name);
//  }
  *(field->fresh_cpu) = YES;
#endif
  return;
}

void Input2DInt_GPU(FieldInt2D *field, int line, const char *string){
  int status;
#ifdef GPU
  if(!(*(field->fresh_gpu))) {
    status = Host2Dev2DInt(field);
    //    printf("%d\n",status);
    //    printf("Copying %s from Dev to Host\n",field->name);
    if(status!=0) {
      printf("Error %d in Host2Dev2DInt()! Field: %s\n",status,field->name);
      printf("called from line %d of file %s\n", line, string);
      exit(EXIT_FAILURE);
    }
  }
//  else {
//    printf("Field %s is up to date on Host.\n",field->name);
//  }
  *(field->fresh_gpu) = YES;
#endif
  return;
}

void Input_CPU(Field *field, int line, const char *string){
  int status, i;
  boolean problem = NO, take_action = NO;
#ifdef GPU
  if (!(*(field->fresh_cpu))) take_action = YES;
  if(take_action) {
    if (*(field->fresh_gpu) == NO) problem = YES;
    if (problem) {
      printf ("Problem on CPU %d: you want to transfer the field %s\n", CPU_Rank, field->name);
      printf ("From device to host in file %s at line %d, but the\n", string, line);
      printf ("field is not fresh everywhere. Here are the values of\n");
      printf ("'fresh' on the mesh, on the inside contour, and outside contour\n");
      INSPECT_INT (*(field->fresh_gpu));
    }
    status = Dev2Host3D(field);
    FullArrayComms++;
    //    printf("%d\n",status);
    //printf("D==>H: %s Level: %d Nx %d Ny %d Pitch_gpu %d\n",field->name,field->level,field->desc->gncell[0],field->desc->gncell[1],field->desc->Pitch_gpu);
    if(status!=0) {
      printf("Error %d in Dev2Host3D()! Field: %s Level: %d\n",status,field->name,Current_Level);
      printf("called from line %d of file %s\n", line, string);
      exit(0);
    }
  }
//  else {
//    printf("Field %s is up to date on Host.\n",field->name);
//  }
  *(field->fresh_cpu) = YES;
#endif
  return;
}

void Input_GPU(Field *field, int line, const char *string){
  int i, status;
  boolean problem=NO, take_action=NO;
#ifdef GPU
  if (!(*(field->fresh_gpu))) take_action = YES;
  if(take_action) {
    if (*(field->fresh_cpu) == NO) problem = YES;
    if (problem) {
      printf ("Problem on CPU %d: you want to transfer the field %s\n", CPU_Rank, field->name);
      printf ("From host to device in file %s at line %d, but the\n", string, line);
      printf ("field is not fresh everywhere. Here are the values of\n");
      printf ("'fresh' on the mesh, on the inside contour, and outside contour\n");
      INSPECT_INT (*(field->fresh_cpu));
    }
    status = Host2Dev3D(field);
    FullArrayComms++;
    //printf("H==>D: %s Level: %d Nx %d Ny %d Pitch_gpu %d\n",field->name,field->level,field->desc->gncell[0],field->desc->gncell[1],field->desc->Pitch_gpu);
    if(status!=0) {
      printf("Error %d in Host2Dev3D()! Field: %s Level: %d\n",status,field->name,Current_Level);
      printf("called from line %d of file %s\n", line, string);
      exit(0);
    }
  }
//  else {
//    printf("Field %s is up to date on Device.\n",field->name);
//  }
  *(field->fresh_gpu) = YES;
#endif
  return;
}

void Output_CPU(Field *field, int line, const char *string){
  #ifdef GPU
  int i;
  field->line_origin = line;
  strncpy(field->file_origin, string, MAXLINELENGTH-1);
  *(field->fresh_cpu) = YES;
  *(field->fresh_gpu) = NO;
  //printf("field %s level %d fresh on CPU\n",field->name,field->level);
  #endif
}

void Output_GPU(Field *field, int line, const char *string){
  #ifdef GPU
  int i;
  field->line_origin = line;
  strncpy(field->file_origin, string, MAXLINELENGTH-1);
  *(field->fresh_cpu) = NO;
  *(field->fresh_gpu) = YES;
  //printf("field %s level %d fresh on GPU\n",field->name,field->level);
  #endif
}


void Output2D_CPU(Field2D *field, int line, const char *string){
  *(field->fresh_cpu) = YES;
  *(field->fresh_gpu) = NO;
}

void Output2D_GPU(Field2D *field, int line, const char *string){
  *(field->fresh_gpu) = YES;
  *(field->fresh_cpu) = NO;
}

void Output2DInt_CPU(FieldInt2D *field, int line, const char *string){
  *(field->fresh_cpu) = YES;
  *(field->fresh_gpu) = NO;
}

void Output2DInt_GPU(FieldInt2D *field, int line, const char *string){
  *(field->fresh_gpu) = YES;
  *(field->fresh_cpu) = NO;
}

void Draft (Field *field, int line, const char *string) {
  strncpy (field->file_origin, string, MAXLINELENGTH-1);
  field->line_origin = line;
}

void WhereIsField(Field *field) {
  int id=0;

  if(*(field->fresh_cpu)){
    printf("Field %s is fresh on the CPU %d\n",field->name,id);
    id=1;
  }
  if(*(field->fresh_gpu))
    printf("Field %s is fresh on the GPU %d\n",field->name,id);
}

void WhereIsFieldInt2D(FieldInt2D *field) {
  int id=0;

  if(*(field->fresh_cpu)){
    printf("Field %s is fresh on the CPU %d\n",field->name,id);
    id=1;
  }
  if(*(field->fresh_gpu))
    printf("Field %s is fresh on the GPU %d\n",field->name,id);
}

void SynchronizeHD () {
#ifdef GPU
  Field *current;
  current = ListOfGrids;
  while (current != NULL) {
    if ((*(current->fresh_cpu) == YES) && (*(current->fresh_gpu) == NO)) {
      Host2Dev3D (current);
      //      printf ("Sending %s to device\n", current->name);
      *(current->fresh_gpu) = YES;
    }
    if ((*(current->fresh_cpu) == NO) && (*(current->fresh_gpu) == YES)) {
      Dev2Host3D (current);
      //printf ("Sending %s to host\n", current->name);
      *(current->fresh_cpu) = YES;
    }
    current = current->next;
  }
#endif
}

void WhereIsWho () {
#ifdef GPU
  Field *current;
  char loc[3];
  current = ListOfGrids;
  loc[1]=0;
  while (current != NULL) {
    if ((*(current->fresh_cpu) == YES) && (*(current->fresh_gpu) == YES))
      loc[0] = 'B';
    if ((*(current->fresh_cpu) == NO) && (*(current->fresh_gpu) == YES))
      loc[0] = 'G';
    if ((*(current->fresh_cpu) == YES) && (*(current->fresh_gpu) == NO))
      loc[0] = 'C';
    if ((*(current->fresh_cpu) == NO) && (*(current->fresh_gpu) == NO))
      loc[0] = '?';
    printf ("%-20s%s\n",current->name, loc);
    current = current->next;
  }
#endif
}

void JSInput_CPU(ScalarField *field, int line, const char *string){
  int status, i;
  boolean problem = NO, take_action = NO;
#ifdef GPU
  if (!(field->fresh_cpu)) take_action = YES;
  if(take_action) {
    if (field->fresh_gpu == NO) problem = YES;
    if (problem) {
      printf ("Problem on CPU %d: you want to transfer the field %s\n", CPU_Rank, field->Name);
      printf ("From device to host in file %s at line %d, but the\n", string, line);
      printf ("field is not fresh everywhere. Here are the values of\n");
      printf ("'fresh' on the mesh, on the inside contour, and outside contour\n");
      INSPECT_INT (field->fresh_gpu);
    }
    status = JSDev2Host3D(field);
    FullArrayComms++;
    //    printf("%d\n",status);
    //printf("D==>H: %s Level: %d\n",field->Name,field->level);
    if(status!=0) {
      printf("Error %d in Dev2Host3D()! Field: %s\n",status,field->Name);
      printf("called from line %d of file %s\n", line, string);
      exit(0);
    }
  }
//  else {
//    printf("Field %s is up to date on Host.\n",field->name);
//  }
  field->fresh_cpu = YES;
#endif
  return;
}

void JSInput_GPU(ScalarField *field, int line, const char *string){
  int i, status;
  boolean problem=NO, take_action=NO;
  /* Below we access directly a nested-mesh structure (ScalarField
*). Its fields fresh_cpu and fresh_gpu are directly booleans, not
*pointers. */
#ifdef GPU
  if (!(field->fresh_gpu)) take_action = YES;
  if(take_action) {
    if (field->fresh_cpu == NO) problem = YES;
    if (problem) {
      printf ("Problem on CPU %d: you want to transfer the field %s\n", CPU_Rank, field->Name);
      printf ("From host to device in file %s at line %d, but the\n", string, line);
      printf ("field is not fresh everywhere. Here are the values of\n");
      printf ("'fresh' on the mesh, on the inside contour, and outside contour\n");
      INSPECT_INT (field->fresh_cpu);
    }
    status = JSHost2Dev3D(field);
    FullArrayComms++;
    //printf("H==>D: %s Level: %d\n",field->Name,field->level);
    if(status!=0) {
      printf("Error %d in Host2Dev3D()! Field: %s\n",status,field->Name);
      printf("called from line %d of file %s\n", line, string);
      exit(0);
    }
  }
//  else {
//    printf("Field %s is up to date on Device.\n",field->name);
//  }
  field->fresh_gpu = YES;
#endif
  return;
}

void JSOutput_CPU(ScalarField *field, int line, const char *string){
  field->fresh_cpu = YES;
  field->fresh_gpu = NO;
  // printf("field %s level %d fresh on CPU\n",field->Name,field->level);
}

void JSOutput_GPU(ScalarField *field, int line, const char *string){
  field->fresh_cpu = NO;
  field->fresh_gpu = YES;
  // printf("field %s level %d fresh on GPU\n",field->Name,field->level);
}

void JVInput_CPU(VectorField *field, int di, int line, const char *string){
  int status, i;
  boolean problem = NO, take_action = NO;
#ifdef GPU
  if (!(field->fresh_cpu[di])) take_action = YES;
  if(take_action) {
    if (field->fresh_gpu[di] == NO) problem = YES;
    if (problem) {
      printf ("Problem on CPU %d: you want to transfer the component %d of vector field %s\n", CPU_Rank, di, field->Name);
      printf ("From device to host in file %s at line %d, but the\n", string, line);
      printf ("field is not fresh everywhere. Here are the values of\n");
      printf ("'fresh' on the mesh, on the inside contour, and outside contour\n");
      INSPECT_INT (field->fresh_gpu[di]);
    }
    status = JVDev2Host3D(field,di);
    FullArrayComms++;
    //    printf("%d\n",status);
    //printf("D==>H: %s[%d] Level: %d\n",field->Name,di,field->level);
    if(status!=0) {
      printf("Error %d in Dev2Host3D()! Field: %s[%d]\n",status,field->Name,di);
      printf("called from line %d of file %s\n", line, string);
      exit(0);
    }
  }
//  else {
//    printf("Field %s is up to date on Host.\n",field->name);
//  }
  field->fresh_cpu[di] = YES;
#endif
  return;
}

 void JVInput_GPU(VectorField *field, int di, int line, const char *string){
  int i, status;
  boolean problem=NO, take_action=NO;
#ifdef GPU
  if (!(field->fresh_gpu[di])) take_action = YES;
  if(take_action) {
    if (field->fresh_cpu[di] == NO) problem = YES;
    if (problem) {
      printf ("Problem on CPU %d: you want to transfer the component %d of vector field %s\n", CPU_Rank, di, field->Name);
      printf ("From host to device in file %s at line %d, but the\n", string, line);
      printf ("field is not fresh everywhere. Here are the values of\n");
      printf ("'fresh' on the mesh, on the inside contour, and outside contour\n");
      INSPECT_INT (field->fresh_cpu[di]);
    }
    status = JVHost2Dev3D(field,di);
    FullArrayComms++;
    //printf("H==>D: %s[%d] Level: %d\n",field->Name,di,field->level);
    if(status!=0) {
      printf("Error %d in Host2Dev3D()!Field: %s[%d]\n",status,field->Name,di);
      printf("called from line %d of file %s\n", line, string);
      exit(0);
    }
  }
//  else {
//    printf("Field %s is up to date on Device.\n",field->name);
//  }
  field->fresh_gpu[di] = YES;
#endif
  return;
}

 void JVOutput_CPU(VectorField *field, int di, int line, const char *string){
  field->fresh_cpu[di] = YES;
  field->fresh_gpu[di] = NO;
  //printf("field %s level %d fresh on CPU\n",field->Name,field->level);
}

 void JVOutput_GPU(VectorField *field, int di, int line, const char *string){
  field->fresh_cpu[di] = NO;
  field->fresh_gpu[di] = YES;
  //printf("field %s level %d fresh on GPU\n",field->Name,field->level);
}
