struct icode {
  char codename[128];
  long index;
};

typedef struct icode Icode;

#ifdef INITCODE
Icode InitLib[512];
#else
extern Icode InitLib[512];
#endif
