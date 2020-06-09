
#ifndef INC__SEM_H_
#define INC__SEM_H_

#include <stdint.h>
#include <pthread.h>
#include <time.h>
#ifdef __APPLE__
#include <dispatch/dispatch.h>
#else
#include <semaphore.h>
#endif

typedef struct osindep_sem_t {
#ifdef __APPLE__
    dispatch_semaphore_t    sem;
    volatile int            val;
    pthread_mutex_t         mtx;
#else
    sem_t                   sem;
#endif
} osindep_sem_t;

static inline int osindep_sem_init(osindep_sem_t *s, int pshared, uint32_t value){
#ifdef __APPLE__
    dispatch_semaphore_t *sem = &s->sem;
    s->val = (int)value;
    pthread_mutex_init(&s->mtx, NULL);
    *sem = dispatch_semaphore_create(value);
    return 0;
#else
    return sem_init(&s->sem, pshared, value);
#endif
}

static inline int osindep_sem_wait(osindep_sem_t *s){
#ifdef __APPLE__
    int ret;
    ret = (int)dispatch_semaphore_wait(s->sem, DISPATCH_TIME_FOREVER);
    if (ret) return ret;
    pthread_mutex_lock(&s->mtx);
    s->val--;
    pthread_mutex_unlock(&s->mtx);
    return 0;
#else
    /*int r;
    do {
        r = sem_wait(&s->sem);
    } while (r == -1 && errno == EINTR);*/
    return sem_wait(&s->sem);
#endif
}

static inline int osindep_sem_trywait(osindep_sem_t *s){
#ifdef __APPLE__
    int ret;
    ret = (int)dispatch_semaphore_wait(s->sem, DISPATCH_TIME_NOW);
    if (ret) return ret;
    pthread_mutex_lock(&s->mtx);
    s->val--;
    pthread_mutex_unlock(&s->mtx);
    return 0;
#else
    /*int r;
    do {
        r = sem_wait(&s->sem);
    } while (r == -1 && errno == EINTR);*/
    /*
    int r = sem_trywait(&s->sem);
    if (r == 0) {
    	std::clog << "TryWait Successful Decrement" << std::endl;
    }
    else {
    	std::clog << "TryWait Fail Decrement" << std::endl;
    }
    return r;*/
    return sem_trywait(&s->sem);
#endif
}

static inline int osindep_sem_waitforduration(osindep_sem_t *s, int64_t ns){
#ifdef __APPLE__
    dispatch_semaphore_wait(s->sem, dispatch_time(DISPATCH_TIME_NOW, ns));
    pthread_mutex_lock(&s->mtx);
    s->val--;
    pthread_mutex_unlock(&s->mtx);
    return 0;
#else
    /*int r;
    do {
        r = sem_wait(&s->sem);
    } while (r == -1 && errno == EINTR);*/
    struct timespec timeout;
    int ret;
    ret = clock_gettime(CLOCK_REALTIME, &timeout);
    if (ret) return ret;
    timeout.tv_nsec += ns % 1000000000;
    timeout.tv_sec += timeout.tv_nsec / 1000000000;
    timeout.tv_sec += ns / 1000000000;
    timeout.tv_nsec %= 1000000000;
    return sem_timedwait(&s->sem, &timeout);
#endif
}

static inline int osindep_sem_post(osindep_sem_t *s){
#ifdef __APPLE__
    pthread_mutex_lock(&s->mtx);
    s->val++;
    pthread_mutex_unlock(&s->mtx);
    dispatch_semaphore_signal(s->sem);
    return 0;
#else
    //std::clog << "SemPost" << std::endl;
    return sem_post(&s->sem);
#endif
}

static inline int osindep_sem_getvalue(osindep_sem_t *s, int *value){
#ifdef __APPLE__
    *value =  s->val;
    return 0;
#else
    return sem_getvalue(&s->sem, value);
#endif
}

static inline int osindep_sem_destroy(osindep_sem_t *s){
#ifdef __APPLE__
    return 0;
#else
    return sem_destroy(&s->sem);
#endif
}

#endif
