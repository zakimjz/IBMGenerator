//
// command line help  & get arguments
//
#include "glob.h"
#include "gen.h"
#include <string.h> //added by MJZaki for strcmp

#define VERSION "Version dated June 3, 1995"

extern char data_file[];
extern char pat_file[];
extern char tax_file[];
//added by MJZaki
extern char ntpc_file[];
extern char conf_file[];

static BOOLEAN userfile = FALSE;

void err_msg(char *str)
{
  cerr << str << flush;
  exit(1);
}


void print_version(void)
{
  cerr << VERSION << endl;
}


void command_line(TaxPar &par)
{
  par.calc_values();

  cerr << "Command Line Options:\n";
  cerr << "  -ntrans number_of_transactions (in 1000's) (default: "
       << par.ntrans/1000 << ")\n";
  cerr << "  -tlen avg_items_per_transaction (default: " << par.tlen << ")\n";
  cerr << "  -nitems number_of_different_items (in '000s) (default: " << par.nitems
       << ")\n";
  cerr << "  -nroots number_of_roots (default: " << par.nroots << ")\n";
  cerr << "  -nlevels number_of_different_levels (default: " << par.nlevels
       << ")\n";
  cerr << "  -fanout average_fanout (default: " << par.fanout << ")\n";
  cerr << "  -depth affects_average_depth_of_items_in_itemsets (default: "
    << par.depth_ratio << ")\n";
  cerr << endl;

  cerr << "  -npats number_of_patterns (default: " << par.lits.npats << ")\n";
  cerr << "  -patlen avg_length_of_maximal_pattern (default: "
       << par.lits.patlen << ")\n";
  cerr << "  -corr correlation_between_patterns (default: " << par.lits.corr
       << ")\n";
  cerr << "  -conf avg_confidence_in_a_rule (default: " << par.lits.conf
       << ")\n";
  cerr << endl;

  cerr << "  -fname <filename> (write to filename.data and filename.pat)\n";
  cerr << "  -ascii (Write data in ASCII format; default: " << (par.ascii? "True": "False") << ")\n";
  cerr << "  -randseed # (reset seed used generate to x-acts; must be negative)\n";
  cerr << "  -version (to print out version info)\n";
  exit(1);
}


void command_line(TransPar &par)
{
  cerr << "Command Line Options:\n";
  cerr << "  -ntrans number_of_transactions (in 1000's) (default: "
       << par.ntrans/1000 << ")\n";
  cerr << "  -tlen avg_items_per_transaction (default: " << par.tlen << ")\n";
  cerr << "  -nitems number_of_different_items (in '000s) (default: " << par.nitems
       << ")\n";
  cerr << endl;

  cerr << "  -npats number_of_patterns (default: " << par.lits.npats << ")\n";
  cerr << "  -patlen avg_length_of_maximal_pattern (default: "
       << par.lits.patlen << ")\n";
  cerr << "  -corr correlation_between_patterns (default: " << par.lits.corr
       << ")\n";
  cerr << "  -conf avg_confidence_in_a_rule (default: " << par.lits.conf
       << ")\n";
  cerr << endl;

  cerr << "  -fname <filename> (write to filename.data and filename.pat)\n";
  cerr << "  -ascii (default: " << (par.ascii? "True": "False") << ")\n";
  cerr << "  -randseed # (reset seed used generate to x-acts; must be negative)\n";
  cerr << "  -version (to print out version info)\n";
  exit(1);
}


void command_line(SeqPar &par)
{
  cerr << "Command Line Options:\n";
  cerr << "  -ncust number_of_customers (in 1000's) (default: "
       << par.ncust/1000 << ")\n";
  cerr << "  -slen avg_trans_per_customer (default: " << par.slen << ")\n";
  cerr << "  -tlen avg_items_per_transaction (default: " << par.tlen << ")\n";
  cerr << "  -nitems number_of_different_items (in '000s) (default: "
       << par.nitems << ")\n";
  cerr << "  -rept repetition-level (default: " << par.rept << ")\n";
  cerr << endl;

  cerr << "  -seq.npats number_of_seq_patterns (default: " << par.lseq.npats
       << ")\n";
  cerr << "  -seq.patlen avg_length_of_maximal_pattern (default: "
       << par.lseq.patlen << ")\n";
  cerr << "  -seq.corr correlation_between_patterns (default: "
       << par.lseq.corr << ")\n";
  cerr << "  -seq.conf avg_confidence_in_a_rule (default: " << par.lseq.conf
       << ")\n";
  cerr << endl;

  cerr << "  -lit.npats number_of_patterns (default: " << par.lits.npats
       << ")\n";
  cerr << "  -lit.patlen avg_length_of_maximal_pattern (default: "
       << par.lits.patlen << ")\n";
  cerr << "  -lit.corr correlation_between_patterns (default: "
       << par.lits.corr << ")\n";
  cerr << "  -lit.conf avg_confidence_in_a_rule (default: " << par.lits.conf
       << ")\n";
  cerr << endl;

  cerr << "  -fname <filename> (write to filename.data and filename.pat)\n";
  cerr << "  -ascii (Write data in ASCII format; default: " << (par.ascii? "True": "False") << ")\n";
  cerr << "  -version (to print out version info)\n";
  exit(1);
}


void cat_fname(char *str1, char *str2)
{
  if (userfile) return;
  
  strcat(data_file, str1);
  strcat(pat_file, str1);
  strcat(tax_file, str1);
  strcat(data_file, str2);
  strcat(pat_file, str2);
  strcat(tax_file, str2);
}


void get_args(TransPar &par, int argc, char **argv)
{
  LINT arg_pos = 2;
  
  strcpy(data_file, "data");
  strcpy(pat_file, "pat");
  strcpy(tax_file, "tax");
  while (arg_pos < argc)
    {
      if (strcmp(argv[arg_pos], "-ntrans") == 0) {
	par.ntrans = (int) (1000 * atof(argv[++arg_pos]));
	cat_fname(".ntrans_", argv[arg_pos]);
	arg_pos++;
	if (par.ntrans < 1) err_msg("ntrans must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-tlen") == 0) {
	par.tlen = atof(argv[++arg_pos]);
	cat_fname(".tlen_", argv[arg_pos]);
	arg_pos++;
	if (par.tlen < 1) err_msg("tlen must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-nitems") == 0) {
	par.nitems = (int) (1000 * atof(argv[++arg_pos]));
	cat_fname(".nitems_", argv[arg_pos]);
	arg_pos++;
	if (par.nitems < 1) err_msg("nitems must be >= 1\n");
	continue;
      }

      else if (strcmp(argv[arg_pos], "-npats") == 0) {
	par.lits.npats = atoi(argv[++arg_pos]);
	cat_fname(".npats_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.npats < 1) err_msg("npats must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-patlen") == 0) {
	par.lits.patlen = atof(argv[++arg_pos]);
	cat_fname(".patlen_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.patlen <= 0) err_msg("patlen must be > 0\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-corr") == 0) {
	par.lits.corr = atof(argv[++arg_pos]);
	cat_fname(".corr_", argv[arg_pos]);
	arg_pos++;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-conf") == 0) {
	par.lits.conf = atof(argv[++arg_pos]);
	cat_fname(".conf_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.conf > 1 || par.lits.conf < 0) 
	  err_msg("conf must be between 0 and 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-fname") == 0) {
        strcpy(data_file, argv[++arg_pos]);
        strcat(data_file, ".data");
        strcpy(pat_file, argv[arg_pos]);
        strcat(pat_file, ".pat");
        //added by MJZaki 
        strcpy(conf_file, argv[arg_pos]);
        strcat(conf_file, ".conf");
        strcpy(ntpc_file, argv[arg_pos]);
        strcat(ntpc_file, ".ntpc");
        arg_pos++;
        userfile = TRUE;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-ascii") == 0) {
	par.ascii = TRUE;
	cat_fname(".ascii", "");
	arg_pos++;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-randseed") == 0) {
	par.seed = atoi(argv[++arg_pos]);
	arg_pos++;
	if (par.seed >= 0)
	  err_msg("randseed must be negative.\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-version") == 0) {
	cout << VERSION << endl;
	exit(0);
      }
      else {
	command_line(par);
      }
    }  // end while
}


void get_args(TaxPar &par, int argc, char **argv)
{
  LINT arg_pos = 2;
  
  strcpy(data_file, "data");
  strcpy(pat_file, "pat");
  strcpy(tax_file, "tax");
  while (arg_pos < argc)
    {
      if (strcmp(argv[arg_pos], "-ntrans") == 0) {
	par.ntrans = (int) (1000 * atof(argv[++arg_pos]));
	cat_fname(".ntrans_", argv[arg_pos]);
	arg_pos++;
	if (par.ntrans < 1) err_msg("ntrans must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-tlen") == 0) {
	par.tlen = atof(argv[++arg_pos]);
	cat_fname(".tlen_", argv[arg_pos]);
	arg_pos++;
	if (par.tlen < 1) err_msg("tlen must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-nitems") == 0) {
	par.nitems = (int) (1000 * atof(argv[++arg_pos]));
	cat_fname(".nitems_", argv[arg_pos]);
	arg_pos++;
	if (par.nitems < 1) err_msg("nitems must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-nroots") == 0) {
	par.nroots = atoi(argv[++arg_pos]);
	cat_fname(".nroots_", argv[arg_pos]);
	arg_pos++;
	if (par.nroots < 1) err_msg("nroots must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-nlevels") == 0) {
	par.nlevels = atof(argv[++arg_pos]);
	cat_fname(".nlevels_", argv[arg_pos]);
	arg_pos++;
	if (par.nlevels < 1) err_msg("nlevels must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-fanout") == 0) {
	par.fanout = atof(argv[++arg_pos]);
	cat_fname(".fanout_", argv[arg_pos]);
	arg_pos++;
	if (par.fanout < 1) err_msg("fanout must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-depth") == 0) {
	par.depth_ratio = atof(argv[++arg_pos]);
	cat_fname(".depth_", argv[arg_pos]);
	arg_pos++;
	if (par.depth_ratio <= 0) err_msg("fanout must be > 0\n");
	continue;
      }

      else if (strcmp(argv[arg_pos], "-npats") == 0) {
	par.lits.npats = atoi(argv[++arg_pos]);
	cat_fname(".npats_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.npats < 1) err_msg("npats must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-patlen") == 0) {
	par.lits.patlen = atof(argv[++arg_pos]);
	cat_fname(".patlen_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.patlen <= 0) err_msg("patlen must be > 0\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-corr") == 0) {
	par.lits.corr = atof(argv[++arg_pos]);
	cat_fname(".corr_", argv[arg_pos]);
	arg_pos++;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-conf") == 0) {
	par.lits.conf = atof(argv[++arg_pos]);
	cat_fname(".conf_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.conf > 1 || par.lits.conf < 0) 
	  err_msg("conf must be between 0 and 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-fname") == 0) {
        strcpy(data_file, argv[++arg_pos]);
        strcat(data_file, ".data");
        strcpy(pat_file, argv[arg_pos]);
        strcat(pat_file, ".pat");
        strcpy(tax_file, argv[arg_pos]);
        strcat(tax_file, ".tax");
        //added by MJZaki 
        strcpy(conf_file, argv[arg_pos]);
        strcat(conf_file, ".conf");
        strcpy(ntpc_file, argv[arg_pos]);
        strcat(ntpc_file, ".ntpc");
        arg_pos++;
        userfile = TRUE;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-ascii") == 0) {
	par.ascii = TRUE;
	cat_fname(".ascii", "");
	arg_pos++;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-randseed") == 0) {
	par.seed = atoi(argv[++arg_pos]);
	arg_pos++;
	if (par.seed >= 0)
	  err_msg("randseed must be negative.\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-version") == 0) {
	cout << VERSION << endl;
	exit(0);
      }
      else {
	command_line(par);
      }
    }  // end while

  par.calc_values();
}


void get_args(SeqPar &par, int argc, char **argv)
{
  LINT arg_pos = 2;
  
  strcpy(data_file, "data");
  strcpy(pat_file, "pat");
  strcpy(tax_file, "tax");
  while (arg_pos < argc)
    {
      if (strcmp(argv[arg_pos], "-ncust") == 0) {
	par.ncust = (int) (1000 * atof(argv[++arg_pos]));
	cat_fname(".ncust_", argv[arg_pos]);
	arg_pos++;
	if (par.ncust < 1) err_msg("ntrans must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-slen") == 0) {
	par.slen = atof(argv[++arg_pos]);
	cat_fname(".slen_", argv[arg_pos]);
	arg_pos++;
	if (par.slen < 1) err_msg("slen must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-tlen") == 0) {
	par.tlen = atof(argv[++arg_pos]);
	cat_fname(".tlen_", argv[arg_pos]);
	arg_pos++;
	if (par.tlen < 1) err_msg("tlen must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-nitems") == 0) {
	par.nitems = (int) (1000 * atof(argv[++arg_pos]));
	cat_fname(".nitems_", argv[arg_pos]);
	arg_pos++;
	if (par.nitems < 1) err_msg("nitems must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-rept") == 0) {
	par.rept = atof(argv[++arg_pos]);
	cat_fname(".rept_", argv[arg_pos]);
	arg_pos++;
	if (par.rept < 0 || par.rept > 1) 
	  err_msg("repetition-level must be between 0 and 1\n");
	continue;
      }

      else if (strcmp(argv[arg_pos], "-seq.npats") == 0) {
	par.lseq.npats = atoi(argv[++arg_pos]);
	cat_fname(".seq.npats_", argv[arg_pos]);
	arg_pos++;
	if (par.lseq.npats < 1) err_msg("npats must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-seq.patlen") == 0) {
	par.lseq.patlen = atof(argv[++arg_pos]);
	cat_fname(".seq.patlen_", argv[arg_pos]);
	arg_pos++;
	if (par.lseq.patlen <= 0) err_msg("patlen must be > 0\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-seq.corr") == 0) {
	par.lseq.corr = atof(argv[++arg_pos]);
	cat_fname(".seq.corr_", argv[arg_pos]);
	arg_pos++;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-seq.conf") == 0) {
	par.lseq.conf = atof(argv[++arg_pos]);
	cat_fname(".seq.conf_", argv[arg_pos]);
	arg_pos++;
	if (par.lseq.conf > 1 || par.lseq.conf < 0) 
	  err_msg("conf must be between 0 and 1\n");
	continue;
      }

      else if (strcmp(argv[arg_pos], "-lit.npats") == 0) {
	par.lits.npats = atoi(argv[++arg_pos]);
	cat_fname(".lit.npats_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.npats < 1) err_msg("npats must be >= 1\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-lit.patlen") == 0) {
	par.lits.patlen = atof(argv[++arg_pos]);
	cat_fname(".lit.patlen_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.patlen <= 0) err_msg("patlen must be > 0\n");
	continue;
      }
      else if (strcmp(argv[arg_pos], "-lit.corr") == 0) {
	par.lits.corr = atof(argv[++arg_pos]);
	cat_fname(".lit.corr_", argv[arg_pos]);
	arg_pos++;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-lit.conf") == 0) {
	par.lits.conf = atof(argv[++arg_pos]);
	cat_fname(".lit.conf_", argv[arg_pos]);
	arg_pos++;
	if (par.lits.conf > 1 || par.lits.conf < 0) 
	  err_msg("conf must be between 0 and 1\n");
	continue;
      }

      else if (strcmp(argv[arg_pos], "-fname") == 0) {
        strcpy(data_file, argv[++arg_pos]);
        strcat(data_file, ".data");
        strcpy(pat_file, argv[arg_pos]);
        strcat(pat_file, ".pat");
        strcpy(tax_file, argv[arg_pos]);
        strcat(tax_file, ".tax");
        //added by MJZaki 
        strcpy(conf_file, argv[arg_pos]);
        strcat(conf_file, ".conf");
        strcpy(ntpc_file, argv[arg_pos]);
        strcat(ntpc_file, ".ntpc");
        arg_pos++;
        userfile = TRUE;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-incr") == 0) {
         par.CustRatio = atof(argv[++arg_pos]);
         par.TransRatio = atof(argv[++arg_pos]);
         par.IncrTid = atoi(argv[++arg_pos]);
         ++arg_pos;
         cout << "INCR :" << par.CustRatio << " " << par.TransRatio 
              << " " << par.IncrTid << endl;
         continue;
      }
      else if (strcmp(argv[arg_pos], "-ascii") == 0) {
	par.ascii = TRUE;
	cat_fname(".ascii", "");
	arg_pos++;
	continue;
      }
      else if (strcmp(argv[arg_pos], "-version") == 0) {
	cout << VERSION << endl;
	exit(0);
      }
      else {
	command_line(par);
      }
    }  // end while
}
