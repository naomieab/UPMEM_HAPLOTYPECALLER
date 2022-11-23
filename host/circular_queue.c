#include "circular_queue.h"

#include <semaphore.h>
#include <pthread.h>
#include <stdbool.h>
#include <malloc.h>
#include <assert.h>

void queue_init(struct queue_t* queue, int size) {
	queue->size = size;
	queue->available = malloc(sizeof(bool)*(size));
	queue->used = malloc(sizeof(bool)*(size));
	for (int i=0; i<size; i++) {
		queue->available[i] = false;
		queue->used[i] = false;
	}
	sem_init(&(queue->used_semaphore), 0, 0);
	sem_init(&(queue->free_semaphore), 0, size-1);
	pthread_mutex_init(&(queue->mutex), NULL);
	queue->next_to_put = 0;
	queue->next_to_release = 0;
	queue->next_to_make_available = 0;
	queue->next_to_take = 0;
    queue->queue_closed = false;
    queue->number_of_writers = 0;
}

void queue_close(struct queue_t* queue, int max_consumers) {
    queue->queue_closed = true;
    for (int i=0; i<max_consumers; i++) {
        sem_post(&queue->used_semaphore);
    }
}

void queue_open_writer(struct queue_t* queue) {
    pthread_mutex_lock(&queue->mutex);
    queue->number_of_writers++;
    pthread_mutex_unlock(&queue->mutex);
}

void queue_close_writer(struct queue_t* queue, int max_consumers) {
    pthread_mutex_lock(&queue->mutex);
    if (--(queue->number_of_writers) <= 0) {
        queue_close(queue, max_consumers);
    }
    pthread_mutex_unlock(&queue->mutex);
}

/*
Ask for an element in the queue (blocking)
Returns the index of the element in the queue that can be used.
*/
int queue_take(struct queue_t* queue) {
	sem_wait(&queue->used_semaphore);
	pthread_mutex_lock(&queue->mutex);
    if (queue->queue_closed && queue->available[queue->next_to_take] == false) {
        pthread_mutex_unlock(&queue->mutex);
        return -1;
    }
	int take_id = queue->next_to_take++;
	queue->next_to_take %= queue->size;
	assert(queue->available[take_id] == true);
	queue->available[take_id] = false;
	pthread_mutex_unlock(&queue->mutex);
	return take_id;
}

/*
Mark an element in the queue as ready to be overwritten.
*/
void queue_release(struct queue_t* queue, int release_id) {
	pthread_mutex_lock(&queue->mutex);
	queue->used[release_id] = false;
	if ((queue->next_to_release)%queue->size == release_id) {
		while (queue->used[queue->next_to_release] == false &&
			   queue->next_to_release != queue->next_to_take) {
			queue->next_to_release+=1;
			queue->next_to_release%=queue->size;
			assert(sem_post(&queue->free_semaphore) == 0);
		}
	}
	pthread_mutex_unlock(&queue->mutex);
}

/*
Ask for ownership of a place in the queue
*/
int queue_put(struct queue_t* queue) {
	sem_wait(&queue->free_semaphore);
	pthread_mutex_lock(&queue->mutex);
	int put_id = queue->next_to_put++;
	queue->next_to_put%=queue->size;
	assert(queue->used[put_id] == false);
	queue->used[put_id] = true;
	pthread_mutex_unlock(&queue->mutex);
	return put_id;
}

/*
Mark the place in the queue as available for taking
(Make sure you own it first)
*/
void queue_make_available(struct queue_t* queue, int available_id) {
	pthread_mutex_lock(&queue->mutex);
	queue->available[available_id] = true;
	if ((queue->next_to_make_available)%queue->size == available_id) {
		while (queue->available[queue->next_to_make_available] &&
			   queue->next_to_make_available != queue->next_to_put) {
			queue->next_to_make_available+=1;
			queue->next_to_make_available%=queue->size;
			assert(sem_post(&queue->used_semaphore) == 0);
		}
	}
	pthread_mutex_unlock(&queue->mutex);
}
