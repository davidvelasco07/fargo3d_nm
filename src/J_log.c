#include "fargo3d.h"
#include "J_jupiter.h"

#define NBMAXMSG 100

static boolean openinfo=NO;
static boolean openwarning=NO;
static FILE *infohdl;
static FILE *warninghdl;
static FILE *errorhdl;
static boolean openerror=NO;
static FILE *errorhdl;
static char info_buf[NBMAXMSG][MAXLINELENGTH];
static char warning_buf[NBMAXMSG][MAXLINELENGTH];
static char error_buf[NBMAXMSG][MAXLINELENGTH];
static long count_info=0;
static long count_warning=0;
static long count_error=0;

void pInfo (const char *template, ...)
{
  va_list ap;
  long i;
  if (CPU_Rank) return;		
/* All that junk below because first infos might be issued while
   'OUTPUTDIR' is not known yet, so that the info file cannot be
   created.  Its creation is allowed (and previous infos are
   flushed) by setting AllowFlushLog to YES */
  if (count_info < NBMAXMSG) {
    va_start (ap, template);
    vsprintf (info_buf[count_info], template, ap);
    va_end (ap);
    count_info++;
  }
  if (AllowFlushLog == YES) {
    if (openinfo == NO) {
      infohdl = prs_open ("info.log");
      openinfo = YES;
    }
    for (i= 0; i < count_info; i++) 
      fprintf (infohdl, "%s", info_buf[i]);
    count_info = 0;
  }
  fflush (infohdl);
}

void pWarning (const char *template, ...)
{
  va_list ap;
  long i;
  if (CPU_Rank) return;		
/* All that  junk below because  first warnings might be  issued while
   'OUTPUTDIR' is  not known yet, so  that the warning  file cannot be
   created.  Its  creation  is  allowed  (and  previous  warnings  are
   flushed) by setting AllowFlushLog to YES */
  if (count_warning < NBMAXMSG) {
    va_start (ap, template);
    vsprintf (warning_buf[count_warning], template, ap);
    va_end (ap);
    count_warning++;
  }
  if (AllowFlushLog == YES) {
    if (openwarning == NO) {
      warninghdl = prs_open ("warning.log");
      openwarning = YES;
      fprintf (stderr, "\n\n #    #    ##    #####   #    #     #    #    #   ####\n");
      fprintf (stderr, " #    #   #  #   #    #  ##   #     #    ##   #  #    #\n");
      fprintf (stderr, " #    #  #    #  #    #  # #  #     #    # #  #  #\n");
      fprintf (stderr, " # ## #  ######  #####   #  # #     #    #  # #  #  ###\n");
      fprintf (stderr, " ##  ##  #    #  #   #   #   ##     #    #   ##  #    #\n");
      fprintf (stderr, " #    #  #    #  #    #  #    #     #    #    #   ####\n");
      fprintf (stderr, "\n\nWARNING file not empty.\n");
      fprintf (stderr, "It is a good idea to check %swarning.log.\n\n", OUTPUTDIR);
    }
    for (i= 0; i < count_warning; i++) 
      fprintf (warninghdl, "%s\n", warning_buf[i]);
    count_warning = 0;
    fflush (warninghdl);
  }
}

//void pError (const char *template, ...)
//{
//  va_list ap;
//  long i;
//  if (CPU_Rank) return;		
///* All that junk below because first errors might be issued while
//   'OUTPUTDIR' is not known yet, so that the error file cannot be
//   created.  Its creation is allowed (and previous errors are
//   flushed) by setting AllowFlushLog to YES */
//  if (count_error < NBMAXMSG) {
//    va_start (ap, template);
//    vsprintf (error_buf[count_error], template, ap);
//    va_end (ap);
//    count_error++;
//  }
//  if (AllowFlushLog == YES) {
//    if (openerror == NO) {
//      errorhdl = prs_open ("error.log");
//      openerror = YES;
//    }
//    for (i= 0; i < count_error; i++) 
//      fprintf (errorhdl, "%s", error_buf[i]);
//    count_error = 0;
//  }
//  fflush (errorhdl);
//}
//
//void FlushLog () {
//  AllowFlushLog = YES;
//  if (count_warning) pWarning ("");
//  if (count_info) pInfo ("");
//  fflush (infohdl);
//  fflush (warninghdl);
//}

//void FlushInfo () {
//  pInfo ("");
//  fflush (infohdl);
//}
