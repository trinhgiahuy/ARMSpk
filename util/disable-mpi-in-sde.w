// fnall iterates over all functions minus the named functions.
// {{fnall <iterator variable name> <function A> <function b> ... }}
// use 'foo' for the required iterator name
{{fnall foo}} {
  __SSC_MARK(0x222);
  {{callfn}}
  //printf("Hello from {{foo}}!\n");
  __SSC_MARK(0x111);
}
{{endfnall}}
