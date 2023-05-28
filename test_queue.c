#include "host/circular_queue.h"

#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#define QUEUE_SIZE 100
#define DATA_SIZE 1000
#define CONSUMER_THREADS 4
#define PRODUCER_THREADS 4

static struct queue_t queue;

int data[QUEUE_SIZE];
int result[DATA_SIZE];

void* producer(void* arg) {
	int start = ((int*) arg)[0];
	int end = ((int*) arg)[1];
	for (int i=start; i<end; i++) {
		int index = queue_put(&queue);
		data[index] = i;
		queue_make_available(&queue, index);
	}
	return 0;
}

void* consumer(void* arg) {
	int n = ((int*) arg)[0];
	for (int i=0; i<n; i++) {
		int index = queue_take(&queue);
        if (index < 0) {
            return 0;
        }
		int result_index = data[index];
		result[result_index] += 1;
		queue_release(&queue, index);
	}
	return 0;
}

int thread_starts[PRODUCER_THREADS+1];
int data_per_consumer = DATA_SIZE/CONSUMER_THREADS+2;

int main() {
	pthread_t threads[8];
	for (int i=0; i<DATA_SIZE; i++) {
		result[i] = 0;
	}
	queue_init(&queue, QUEUE_SIZE);
	for (int i=0; i<PRODUCER_THREADS; i++) {
		thread_starts[i] = DATA_SIZE/PRODUCER_THREADS*i;
	}
	thread_starts[PRODUCER_THREADS] = DATA_SIZE;
	for (int i=0; i<CONSUMER_THREADS; i++) {
		pthread_create(&threads[PRODUCER_THREADS+i], NULL, consumer, &data_per_consumer);
	}
	for (int i=0; i<PRODUCER_THREADS; i++) {
		pthread_create(&threads[i], NULL, producer, &thread_starts[i]);
	}
	for (int i=0; i<PRODUCER_THREADS; i++) {
		pthread_join(threads[i], NULL);
	}
    queue_close(&queue, CONSUMER_THREADS);
	for (int i=PRODUCER_THREADS; i<PRODUCER_THREADS+CONSUMER_THREADS; i++) {
		pthread_join(threads[i], NULL);
	}
	bool error = false;
	for (int i=0; i<DATA_SIZE; i++) {
		if (result[i] != 1) {
			fprintf(stderr, "error at index: %d (%d mappings)\n", i, result[i]);
			error = true;
		}
	}
	if (error) {
		return 1;
	}
	fprintf(stderr, "all good\n");
	return 0;
}
