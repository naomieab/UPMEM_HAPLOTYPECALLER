#ifndef __LOG_H__
#define __LOG_H__

#include <stdio.h>

#define LOG_LVL_NOT    0
#define LOG_LVL_FATAL  1
#define LOG_LVL_ERROR  2
#define LOG_LVL_WARN   3
#define LOG_LVL_INFO   4
#define LOG_LVL_DEBUG  5
#define LOG_LVL_TRACE  6

#define LOG_LEVEL LOG_LVL_INFO
#define LOG_FILE  stderr


#define LOG(lvl, ...) LOG_ ## lvl ( ## __VA_ARGS)

#if LOG_LEVEL >= LOG_LVL_FATAL
#define LOG_FATAL(...) fprintf(LOG_FILE, ## __VA_ARGS__)
#else
#define LOG_FATAL(...)
#endif

#if LOG_LEVEL >= LOG_LVL_ERROR
#define LOG_ERROR(...) fprintf(LOG_FILE, ## __VA_ARGS__)
#else
#define LOG_ERROR(...)
#endif

#if LOG_LEVEL >= LOG_LVL_WARN
#define LOG_WARN(...) fprintf(LOG_FILE, ## __VA_ARGS__)
#else
#define LOG_WARN(...)
#endif

#if LOG_LEVEL >= LOG_LVL_INFO
#define LOG_INFO(...) fprintf(LOG_FILE, ## __VA_ARGS__)
#else
#define LOG_INFO(...)
#endif

#if LOG_LEVEL >= LOG_LVL_DEBUG
#define LOG_DEBUG(...) fprintf(LOG_FILE, ## __VA_ARGS__)
#else
#define LOG_DEBUG(...)
#endif

#if LOG_LEVEL >= LOG_LVL_TRACE
#define LOG_TRACE(...) fprintf(LOG_FILE, ## __VA_ARGS__)
#else
#define LOG_TRACE(...)
#endif



#endif //__LOG_H__
