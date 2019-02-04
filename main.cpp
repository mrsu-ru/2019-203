#include <cstdlib>
#include <cstring>
#include "lab.h"
#include <iostream>
#include "alekseev.h"
#include "borisovaem.h"
#include "viktorovichev.h"
#include "volkovda.h"
#include "dryginaea.h"
#include "epifanovats.h"
#include "eremkinnv.h"
#include "igaykiniv.h"
#include "kemaykinais.h"
#include "kiselevis.h"
#include "kobyzevkv.h"
#include "makarovaaa.h"
#include "malovava.h"
#include "maslovma.h"
#include "medvedevama.h"
#include "murtazinrr.h"
#include "nachinkinaoa.h"
#include "rusaevon.h"
#include "sirkinaa.h"
#include "sovetnikovakv.h"
#include "tenishevir.h"
#include "tereshkinki.h"
#include "tokarevaaa.h"
#include "shmelevaov.h"
#include "zhalninrv.h"


void print_usage(char* name);


int main(int argc, char** argv)
{
  if (argc < 3) {
    print_usage(argv[0]);
    return 0;
  }

  lab *l = NULL;
  if (strcmp(argv[1], "alekseev") == 0) {
    l = new alekseev();
  }
  else if (strcmp(argv[1], "borisovaem") == 0) {
    l = new borisovaem();
  }
  else if (strcmp(argv[1],"viktorovichev")==0) {
    l = new viktorovichev();
  }
  else if (strcmp(argv[1],"volkovda")==0) {
    l = new volkovda();
  }
  else if (strcmp(argv[1],"dryginaea")==0) {
    l = new dryginaea();
  }
  else if (strcmp(argv[1],"epifanovats")==0) {
    l = new epifanovats();
  }
  else if (strcmp(argv[1],"eremkinnv")==0) {
    l = new eremkinnv();
  }
  else if (strcmp(argv[1],"igaykiniv")==0) {
    l = new igaykiniv();
  }
  else if (strcmp(argv[1],"kemaykinais")==0) {
    l = new kemaykinais();
  }
  else if (strcmp(argv[1],"kiselevis")==0) {
    l = new kiselevis();
  }
  else if (strcmp(argv[1],"kobyzevkv")==0) {
    l = new kobyzevkv();
  }
  else  if (strcmp(argv[1], "makarovaaa") == 0) {
    l = new makarovaaa();
  }
  else  if (strcmp(argv[1], "malovava") == 0) {
    l = new malovava();
  }
  else  if (strcmp(argv[1], "maslovma") == 0) {
    l = new maslovma();
  }
  else  if (strcmp(argv[1], "medvedevama") == 0) {
    l = new medvedevama();
  }
  else  if (strcmp(argv[1], "murtazinrr") == 0) {
    l = new murtazinrr();
  }
  else  if (strcmp(argv[1], "nachinkinaoa") == 0) {
    l = new nachinkinaoa();
  }
  else  if (strcmp(argv[1], "rusaevon") == 0) {
    l = new rusaevon();
  }
  else  if (strcmp(argv[1], "sirkinaa") == 0) {
    l = new sirkinaa();
  }
  else  if (strcmp(argv[1], "sovetnikovakv") == 0) {
    l = new sovetnikovakv();
  }
  else  if (strcmp(argv[1], "tenishevir") == 0) {
    l = new tenishevir();
  }
  else  if (strcmp(argv[1], "tereshkinki") == 0) {
    l = new tereshkinki();
  }
  else  if (strcmp(argv[1], "tokarevaaa") == 0) {
    l = new tokarevaaa();
  }
  else  if (strcmp(argv[1], "shmelevaov") == 0) {
    l = new shmelevaov();
  }
  else  if (strcmp(argv[1], "zhalninrv") == 0) {
    l = new zhalninrv();
  }
  else  {
    print_usage(argv[0]);
    return 0;
  }

  l->read_file();
  l->run(atoi(argv[2]));
  l->write_result();
  l->check_result();

  //delete l; // TODO:
  return 0;
}


void print_usage(char* name)
{
  std::cout << "Usage:\n\n  " << name << " <fio> <lab_number>\n";
}
