#include <fstream>
#include <stdio.h>
#include "glob.h"
#include "gen.h"
#include <string.h> //added by MJZaki for strcmp

using namespace std;

// prototypes
void command_line(TransPar &par);
void get_args(TransPar &par, int argc, char **argv);
void gen_rules(TransPar &par);
Transaction *mk_tran(StringSetIter &lits, LINT tlen, LINT NPATS,
                     Taxonomy *tax = NULL); //MJZ: added NPATS

void command_line(TaxPar &par);
void get_args(TaxPar &par, int argc, char **argv);
void gen_taxrules(TaxPar &par);

void command_line(SeqPar &par);
void get_args(SeqPar &par, int argc, char **argv);
void print_version(void);

void gen_seq(SeqPar &par);
CustSeq *mk_seq(Cid cid, StringSetIter &lseq, StringSet &lits, LINT slen, LINT tlen);

char data_file[256];
char pat_file[256];
char tax_file[256];
//added by MJZaki
char ntpc_file[256];
char conf_file[256];


void memory_err(void)
{
  cout << "A memory allocation error occurred. \n";
  exit(1);
}


int main(int argc, char **argv)
{

  set_new_handler(memory_err);

  if (strcmp(argv[1], "lit") == 0) {
    // For Rules
    TransPar par;

    get_args(par, argc, argv);   // get arguments
    gen_rules(par);              // generate rules (really, just transactions)
  }

  else if (strcmp(argv[1], "seq") == 0) {
    // For Sequences
    SeqPar par;

    get_args(par, argc, argv);   // get arguments
    gen_seq(par);                // generate sequences
  }

  else if (strcmp(argv[1], "tax") == 0) {
    // For Rules with Taxonomies
    TaxPar par;

    get_args(par, argc, argv);   // get arguments
    gen_taxrules(par);           // generate rules (really, just transactions)
  }

  else if (strcmp(argv[1], "-version") == 0) {
    print_version();
    return 0;
  }

  else {
    cerr << "Synthetic Data Generation, ";
    print_version();
    cerr << "Usage:  " << argv[0] << " lit|tax|seq [options]\n";
    cerr << "        " << argv[0]
      << " lit|tax|seq -help     For more detailed list of options\n";
    return 1;
  }

  return 0;
}


// Generate Transactions
//
void gen_rules(TransPar &par)
{
  StringSet *lits;
  StringSetIter *patterns;
  Transaction *trans;
  PoissonDist *tlen;

  ofstream data_fp;
  ofstream pat_fp;
  ofstream conf_fp; //added by MJZaki

  data_fp.open(data_file);
  pat_fp.open(pat_file);
  conf_fp.open(conf_file); //added by MJZaki
  if (data_fp.fail() || pat_fp.fail() || conf_fp.fail()) {
    cerr << "Error opening output file" << endl;
    exit(1);
  }

  lits = new StringSet(par.nitems, par.lits);

  // Reset random seed generator for before generating transactions
  if (par.seed < 0) RandSeed::set_seed(par.seed);

  tlen = new PoissonDist(par.tlen-1);

  par.write(pat_fp);
  lits->display(pat_fp);

  patterns = new StringSetIter(*lits);
  LINT NTRANS=0;
  //Transaction::set_print_cid(FALSE); // added by me to suppress cid
  for (LINT i = 0; i < par.ntrans; i ++)
    {
      trans = mk_tran(*patterns, (*tlen)()+1,par.lits.npats);
      if (trans->tid < par.mintid) par.mintid = trans->tid;
      if (trans->tid > par.maxtid) par.maxtid = trans->tid;
      if (par.ascii)
	trans->write_asc(data_fp);
      else
	trans->write(data_fp);
      //cout << "TRANS SZ " << trans->size() << endl;
      if (trans->size() > 0) NTRANS++; //added by me: repeat if trans empty
      else i--;
      delete trans;
    }

  data_fp.close();
  pat_fp.close();

  //added by MJZaki
  if (par.ascii){
     conf_fp << NTRANS << "\n";
     conf_fp << par.nitems << "\n";
     conf_fp << par.tlen << "\n";
     conf_fp << par.mintid << "\n";
     conf_fp << par.maxtid << "\n";
  }
  else{
     cout << "WRITING " << NTRANS << " "<< par.nitems << endl;
     conf_fp.write((char *)&NTRANS, sizeof(LINT));
     conf_fp.write((char *)&par.nitems, sizeof(LINT));
     conf_fp.write((char *)&par.tlen, sizeof(FLOAT));
     conf_fp.write((char *)&par.mintid, sizeof(LINT));
     conf_fp.write((char *)&par.maxtid, sizeof(LINT));
  }
  conf_fp.close();
}


// Generate Transactions and Taxonomy
//
void gen_taxrules(TaxPar &par)
{
  Taxonomy *tax;
  StringSet *lits;
  StringSetIter *patterns;
  Transaction *trans;
  PoissonDist *tlen;

  ofstream data_fp;
  ofstream pat_fp;
  ofstream tax_fp;
  ofstream conf_fp; //added by MJZaki

  data_fp.open(data_file);
  pat_fp.open(pat_file);
  tax_fp.open(tax_file);
  conf_fp.open(conf_file); //added by MJZaki

  if (data_fp.fail() || pat_fp.fail() || tax_fp.fail() || conf_fp.fail()) {
    cerr << "Error opening output file" << endl;
    exit(1);
  }

  // generate taxonomy and write it to file
  tax = new Taxonomy(par.nitems, par.nroots, par.fanout, par.depth_ratio);
  if (par.ascii)
    tax->write_asc(tax_fp);
  else
    tax->write(tax_fp);

  tlen = new PoissonDist(par.tlen-1);

  lits = new StringSet(par.nitems, par.lits, tax);

  par.write(pat_fp);
  lits->display(pat_fp);

  patterns = new StringSetIter(*lits);
  LINT NTRANS=0;
  for (LINT i = 0; i < par.ntrans; i ++)
    {
      trans = mk_tran(*patterns, (*tlen)()+1, par.lits.npats, tax);
      if (par.ascii)
	trans->write_asc(data_fp);
      else
	trans->write(data_fp);
      if (trans->size() > 0) NTRANS++;//added by me: repeat if trans empty
      else i--;
      delete trans;
      delete trans;
    }

  data_fp.close();
  pat_fp.close();
  tax_fp.close();

  //added by MJZaki
  if (par.ascii){
     conf_fp << NTRANS << "\n";
     conf_fp << par.nitems << "\n";
     conf_fp << par.tlen << "\n";
  }
  else{
     conf_fp.write((char *)&NTRANS, sizeof(LINT));
     conf_fp.write((char *)&par.nitems, sizeof(LINT));
     int t = (int) par.tlen;
     conf_fp.write((char *)&t, sizeof(LINT));
  }
  conf_fp.close();
}


// Generate a transaction
//
Transaction *mk_tran(StringSetIter &lits,  	// table of patterns
		     LINT tlen,			// transaction length
                     LINT NPATS,
		     Taxonomy *tax
		     )
{
  Transaction *trans;
  StringP pat;

  if (tlen > Transaction::MAXNITEMS)
     tlen = Transaction::MAXNITEMS; //MJZ: can't exceed nitems
  trans = new Transaction(tlen);
  LINT patcnt=0; //MJZ
  while (trans->size() < tlen && patcnt < NPATS)
    {
       patcnt++;
       //cout << trans->size() << " " << tlen << " HERE2\n";
      pat = lits.get_pat();		// get a pattern
      if ( !trans->add(*pat) ) {
	// this pattern didn't fit in the transaction
	lits.unget_pat();
	break;
      }
    }
  return trans;
}


// Generate Sequences
//
void gen_seq(SeqPar &par)
{
  StringSet *lseq;	// potentially large sequences
  StringSetIter *patterns;
  StringSet *lits;	// potentially large itemsets
  CustSeq *cust;	//
  PoissonDist *slen;
  PoissonDist *tlen;

  ofstream data_fp;
  ofstream pat_fp;
  ofstream conf_fp; //added by MJZaki
  ofstream ntpc_fp; //added by MJZaki
  srand48(0);


  data_fp.open(data_file);
  pat_fp.open(pat_file);
  conf_fp.open(conf_file); //added by MJZaki
  ntpc_fp.open(ntpc_file); //added by MJZaki
  LINT *NTPC = new LINT[par.ncust]; //added by MJZaki
  LINT tottrans=0;
  if (data_fp.fail() || pat_fp.fail() || ntpc_fp.fail() || conf_fp.fail()) {
    cerr << "Error opening output file" << endl;
    exit(1);
  }

  slen = new PoissonDist(par.slen-1);
  tlen = new PoissonDist(par.tlen-1);

  lits = new StringSet(par.nitems, par.lits);
  lseq = new StringSet(par.lits.npats, par.lseq, NULL, par.rept, par.rept_var);

//  pat_fp << "Large Itemsets:" << endl;
//  lits->write(pat_fp);
//  pat_fp << endl << endl << "Sequences:" << endl;
  par.write(pat_fp);
  lseq->display(pat_fp, *lits);

  patterns = new StringSetIter(*lseq);
  LINT NCUST=0;
  LINT i;
  for (i = 0; i < par.ncust; i ++)
    {
       if (i%1000 == 0)
          cout << "DONE " << i << endl;
      cust = mk_seq(i+1, *patterns, *lits, (*slen)()+1, (*tlen)()+1);
      if (cust->cid < par.mincustid) par.mincustid = cust->cid;
      if (cust->cid > par.maxcustid) par.maxcustid = cust->cid;

      if (par.ascii)
         NTPC[NCUST] = cust->write_asc(data_fp);
      else
         NTPC[NCUST] = cust->write(data_fp);
      tottrans += NTPC[NCUST];
      if (NTPC[NCUST] > 0) NCUST++;//added by MJZaki: repeat if trans empty
      else i--;
      delete cust;
    }

  data_fp.close();
  pat_fp.close();

  //added by MJZaki
  if (par.ascii){
     conf_fp << NCUST << "\n";
     conf_fp << par.nitems << "\n";
     conf_fp << par.slen << "\n";
     conf_fp << par.tlen << "\n";
     conf_fp << tottrans << "\n";
     conf_fp << par.mincustid << "\n";
     conf_fp << par.maxcustid << "\n";
  }
  else{
     conf_fp.write((char *)&NCUST, sizeof(LINT));
     conf_fp.write((char *)&par.nitems, sizeof(LINT));
     conf_fp.write((char *)&par.slen, sizeof(FLOAT));
     conf_fp.write((char *)&par.tlen, sizeof(FLOAT));
     conf_fp.write((char *)&tottrans, sizeof(LINT));
     conf_fp.write((char *)&par.mincustid, sizeof(LINT));
     conf_fp.write((char *)&par.maxcustid, sizeof(LINT));
  }
  conf_fp.close();

  if (par.ascii){
     for (i=0; i < NCUST; i++)
        ntpc_fp << NTPC[i] << " ";
     ntpc_fp << endl;
  }
  else ntpc_fp.write((char *)NTPC, NCUST*sizeof(LINT));
  ntpc_fp.close();
  delete [] NTPC;
}


// Generate a customer-sequence
//
CustSeq *mk_seq(Cid cid,		// customer-id
		StringSetIter &lseq,	// table of large sequences
		StringSet &lits,	// table of large itemsets
		LINT slen,		// sequence length
		LINT tlen		// avg. transaction length
		)
{
  CustSeq *cust;
  StringP pat;

  cust = new CustSeq(cid, slen, tlen);
  while (cust->size() < slen * tlen)
    {
      pat = lseq.get_pat();      // get a pattern
      if ( !cust->add(*pat, lits) ) {	// transaction full
	lseq.unget_pat();
	break;
      }
    }
  return cust;
}


