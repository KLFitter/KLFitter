#ifndef KLFITTER_UNUSED
#if defined(__GNUC__) 
#define KLFITTER_UNUSED(x) UNUSED_ ## x __attribute__((unused))
#else
#define KLFITTER_UNUSED(x) x
#endif
#endif
