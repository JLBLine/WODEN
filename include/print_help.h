/*! \file print_help.h
  Command line help for the main `woden` executable - this is for when you run
  `woden` directly and only create a binary, rather than using the
  `run_woden.py` script to run `woden` and convert the outputs into `.uvfits`
*/

/**
@brief Print out how to run `woden`, and a run down of all necessary and
optional inputs in the controlling `.json` file

*/
void print_cmdline_help();
