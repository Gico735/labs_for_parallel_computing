#include <stdio.h>
#include <semaphore.h>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

#define SMOKERS 3

sem_t sTobacco, sPaper, sMatch;
sem_t sBartender;

int Tobacco = 0;
int Paper = 0;
int Match = 0;






void *smoker(void * arg)
{
  int id = *(int*)arg;
  for(;;) {
    switch(id) {
      case 0://tobacco
        sem_wait(&sTobacco);
        if(Paper > 0 && Match > 0) {
          Paper = 0;
          Match = 0;
          printf("%d: Got paper and match\n", id);
          sleep(1);
          printf("%d: start smoke\n", id);
          sem_post(&sBartender);
          sleep(5);
          printf("%d:  smoked\n", id);
        } else printf("%d: Error: Paper or match is not avalible!\n", id);
        break;
      case 1://paper
        sem_wait(&sPaper);
        if(Tobacco > 0 && Match > 0) {
          Tobacco = 0;
          Match = 0;
          printf("%d: Got tobacco and match\n", id);
          sleep(1);
          printf("%d: start smoke\n", id);
          sem_post(&sBartender);
          sleep(5);
          printf("%d:  smoked\n", id);
        } else printf("%d: Error: Tobacco or match is not avalible!\n", id);
        break;
      case 2://Match
        sem_wait(&sMatch);
        if(Tobacco > 0 && Paper > 0) {
          Tobacco = 0;
          Paper = 0;
          printf("%d: Got tobacco and paper\n", id);
          sleep(1);
          printf("%d: start smoke\n", id);
          sem_post(&sBartender);
          sleep(5);
          printf("%d:  smoked\n", id);
        } else printf("%d: Error: Tobacco or paper is not avalible!\n", id);
        break;
      default:
        printf("ERROR: Invalid argument ID: %d\n", id);
        exit(1);
        break;
    }
  }
}

void *bartender(void * arg) {
  for(;;) {
    sem_wait(&sBartender); 
    usleep(1000000);  
    int id = rand() % 3;
    switch(id) {
      case 0: 
        puts("Bartender set paper and match");
        Paper++;
        Match++;
        sem_post(&sTobacco);
        break;
      case 1: 
        puts("Bartender set tobacco and match");
        Tobacco++;
        Match++;
        sem_post(&sPaper);
        break;
      case 2: 
        puts("Bartender set tobacco and paper");
        Tobacco++;
        Paper++;
        sem_post(&sMatch);
        break;
      default:
        printf("ERROR wrong Bartender ID: %d\n", id);
        break;
    }
  }
}

int main(int argc, char* argv[]) {
  sem_init(&sTobacco, 0, 0);
  sem_init(&sPaper, 0, 0);
  sem_init(&sMatch, 0, 0);

  sem_init(&sBartender, 0, 1);

  pthread_t tSmokers[SMOKERS], tBartender;

  int smokersID[SMOKERS];

  for(int i = 0; i < SMOKERS; i++) {
    smokersID[i] = i;
    pthread_create(&tSmokers[i], NULL, smoker, &smokersID[i]);
  }
  puts("Smokers created");

  pthread_create(&tBartender, NULL, bartender, NULL);
  puts("Bartender created");

  for (int i = 0; i < SMOKERS; i++)
    pthread_join(tSmokers[i], NULL);

  pthread_join(tBartender, NULL);
  exit(0);
}