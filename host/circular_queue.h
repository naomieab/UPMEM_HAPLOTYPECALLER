#include <semaphore.h>
#include <pthread.h>
#include <stdbool.h>

#ifndef CIRCULAR_QUEUE_H__
#define CIRCULAR_QUEUE_H__

struct queue_t {
	pthread_mutex_t mutex;
	sem_t           used_semaphore;
	sem_t           free_semaphore;
	int             next_to_put;
	int             next_to_release;
	int             next_to_take;
	int             next_to_make_available;
	int             size;
    int             number_of_writers;
	bool*           available;
	bool*           used;
    bool            queue_closed;
};

void queue_init(struct queue_t* queue, int size);
void queue_close(struct queue_t* queue, int max_consumers);
void queue_open_writer(struct queue_t* queue);
void queue_close_writer(struct queue_t* queue, int max_consumers);
int queue_take(struct queue_t* queue);
void queue_release(struct queue_t* queue, int release_id);
int queue_put(struct queue_t* queue);
void queue_make_available(struct queue_t* queue, int available_id);

#endif // CIRCULAR_QUEUE_H__
