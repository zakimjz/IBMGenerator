#include "glob.h"
#include "dist.h"
#include <limits.h>
#include <iostream>
#include <fstream>

using namespace std;


#define DEFNITEMS 10000
//=============================  Parameters  =============================


// Parameters used for StringSet
// StringSet can be either large itemsets, or sequences
//
class PatternPar {
public:
  LINT npats;	// number of patterns
  FLOAT patlen; // average length of pattern
  FLOAT corr;	// correlation between consecutive patterns
  FLOAT conf;	// average confidence in a rule
  FLOAT conf_var;  // variation in the confidence
  LINT seed;

  PatternPar(void)
    : npats(10000), patlen(4.0), corr(0.25), conf(0.75), conf_var(0.1), seed(0)
  {}

  void write(ostream &fp);
};


// Parameters used to generate transactions
//
class TransPar {
public:
  LINT ntrans;	// number of transactions in database
  FLOAT tlen;	// average transaction length
  LINT nitems;	// number of items
  PatternPar lits;	// parameters for potentially large itemsets
  BOOLEAN ascii;        // Generate in ASCII format
  LINT seed;    // Seed to initialize RandSeed with before x-act generation

   LINT mintid; //MJZaki 17/10/02
   LINT maxtid; //MJZaki 17/10/02

   TransPar(void)
      : ntrans(1000000), tlen(10), nitems(100000), ascii(FALSE),
      seed(INIT_SEED), mintid(INT_MAX), maxtid(INT_MIN)
      {}

  void write(ostream &fp);
};


// Parameters used to generate transactions
//
class TaxPar:public TransPar {
public:
  LINT nroots;	 // number of roots
  FLOAT fanout;	 // average fanout at each interiori node
  FLOAT nlevels; // average number of levels
  FLOAT depth_ratio;	 // affects ratio of itemsets chosen from higher levels

  TaxPar(void)
    : nroots(0), fanout(0), nlevels(0), depth_ratio(1)
  {}

  void calc_values(void);	// calculates nroots, given nlevels
  	// default values: nroots = 250, fanout = 5
  void write(ostream &fp);
};

// Parameters used to generate sequences
//
class SeqPar {
public:
  LINT ncust;	// number of customers in database
  FLOAT slen;	// average sequence length
  FLOAT tlen;	// average transaction length
  LINT nitems;	// number of items

  FLOAT rept;		// repetition-level (between 0 and 1)
  FLOAT rept_var;	// variation in repetition-level

  BOOLEAN ascii;        // Generate in ASCII format

  PatternPar lits;	// parameters for potentially large itemsets
  PatternPar lseq;	// parameters for potentially large sequences

   // Static parameters for incremental DB generation
   static FLOAT CustRatio; //what % of customers have increments
   static FLOAT TransRatio; //what % of trans make up the incr (def min=2 trans)
   static LINT IncrTid; //Starting Tid for incr DB

   LINT mincustid; //MJZaki 17/10/02
   LINT maxcustid; //MJZaki 17/10/02

  SeqPar(void)
    : ncust(100000), slen(10), tlen(2.5), nitems(DEFNITEMS),
     rept(0), rept_var(0.1), ascii(FALSE),
     mincustid(INT_MAX), maxcustid(INT_MIN) //MJZaki
  {
    lits.npats = 25000;
    lseq.npats = 5000;
    lits.patlen = 1.25;
    lseq.patlen = 4.0;
  }

  void write(ostream &fp);
};


//------------------------------ Taxonomy ------------------------------


//
// models taxonomy over items as a tree
// 0 is a valid item here (get rid of it when actually adding item
//
class Taxonomy
{
  friend class TaxStat;
private:
  LINT nitems;	// number of items
  LINT nroots;	// number of roots
  FLOAT depth;	// used when assigning probabilities to items

  LINT *par;
  LINT *child_start;
  LINT *child_end;

  static const LINT item_len;  // ASCII field-width of item-id
public:
  Taxonomy(LINT nitems, LINT nroots, FLOAT fanout, FLOAT depth_ratio);
  ~Taxonomy(void);

  void write(ofstream &fp);		// write taxonomy to file
  void write_asc(ofstream &fp);	// write taxonomy to ASCII file
  void display(ofstream &fp);	// display taxonomy (for user)

  FLOAT depth_ratio(void) { return depth; }
  LINT num_roots(void) { return nroots; }
  LINT root(Item itm) { return (par[itm] == -1); }

  LINT num_children(Item itm) { return child_end[itm] -  child_start[itm]; }
  LINT child(Item itm, LINT n) { return child_start[itm]+n; }
	// returns the n'th child of itm
  LINT first_child(Item itm) { return child_start[itm]; }
  LINT last_child(Item itm) { return child_end[itm]-1; }

  Item parent(Item itm) { return par[itm]; }	// -1 => no parent
};


//--------------------------------------------------------------------------


//
// 0 is a valid item here (get rid of it when actually adding item
//
class ItemSet
{
private:
  LINT nitems;		// number of items
  Taxonomy *tax;	// taxonomy (optional)
  FLOAT *cum_prob;	// cumulative probability
  FLOAT *tax_prob;	// cumulative probability of choosing a child
  UniformDist rand;

  void normalize(FLOAT prob[], LINT low, LINT high);
public:
  ItemSet(LINT nitems, Taxonomy *tax = NULL);
  ~ItemSet();
  void display(ofstream &fp);
  Item get_item(void);		// returns a random item (weighted)
  Item specialize(Item itm);	// if no taxonomy, returns itm
  FLOAT weight(Item itm);	// returns prob. of choosing item
};


class StringSet;
class String
{
friend class StringSet;
  LINT nitems;	// number of items
  Item *items;  // list of the items
//  FLOAT *rval;	// random value (used to get random ordering of the items)
//  Item *ritems;	// randomly chosen items
  FLOAT prob;	// probability that this string is chosen
  FLOAT conf;	// probability that this string is corrrupted

//  void shuffle(void);	// shuffles items in string
public:
  String(LINT nitems);
  ~String(void);

  void display(ofstream &fp, LINT prob_comp = 1);
  void display(ofstream &fp, StringSet &lits, LINT prob_comp = 1);
	// prob is multiplied by prob_comp before being writeed

  LINT size(void) { return nitems;}
  Item item(LINT n) { return items[n];} // return nth item of the string
  FLOAT conf_lvl(void) { return conf; }

  void set_size(LINT newsize) { nitems = newsize;}
  void set_item(LINT n, Item itm) { items[n] = itm;}
  void set_conf_lvl(FLOAT newconf) { conf = newconf; }

//  void shuffle(LINT k);	// allows selection of k random items from the string
//  Item rand_item(LINT n) { return ritems[n];} // works with shuffle
};

typedef String *StringP;

class StringSet
{
  friend class StringSetIter;
private:
  ItemSet *items;
  Taxonomy *tax;
  LINT npats;		// number of patterns
  StringP *pat;		// array of patterns
  StringP answer;
  FLOAT *cum_prob;	// cumulative probabilities of patterns

  StringP specialize(LINT i);	// specializes pattern #i
public:
  StringSet(LINT nitems,	  // number of items
	    PatternPar par,	  // npats, patlen, corr, conf & conf_var
	    Taxonomy *tax = NULL, // taxonomy (optional)
	    FLOAT rept = 0,	  // repetition-level
	    FLOAT rept_lvl = 0.2  // variation in repetition-level
	    );
  ~StringSet();
  void display(ofstream &fp);	// for large itemsets
  void display(ofstream &fp, StringSet &lit);	// for sequences
  StringP get_pat(LINT i);	// returns pattern #i
};

class StringSetIter
{
private:
  UniformDist rand;
  StringSet *strset;
  LINT last_pat;	// if -ve, unget_pat() was called
public:
  StringSetIter(StringSet &str_set) : strset(&str_set), last_pat(0) {};
  StringP get_pat(void);	// returns a random pattern
  void unget_pat(void);		// the last pattern is put back in the sequence
};





//--------------------------------------------------------------------------

class Transaction {
private:
  LINT tlen;	// expected number of items in transaction
  LINT nitems;	// number of items currently in transaction
  LINT maxsize;	// size of array items
  LINT *items;	// items in the transaction


  static const LINT cid_len;   // ASCII field-width of customer-id
  static const LINT tid_len;   // ASCII field-width of transaction-id
  static const LINT item_len;  // ASCII field-width of item-id
   static BOOLEAN print_cid;    //MJZ: print cid or not



  void sort(void);
  BOOLEAN add_item(LINT itm);// returns TRUE if added, FALSE if already present
public:
  static LINT tid;	// transaction-id
   static LINT MAXNITEMS;      //MJZ: max number of possible items
  Transaction(LINT sz);
  ~Transaction();

  BOOLEAN add(String &pat, BOOLEAN corrupt = TRUE);
	// adds pattern to transaction
	// returns TRUE if added, FALSE if trans. full
   void write(ofstream &fp, LINT cid = 0, LINT incrtid=-1); //added incrtid --MJZaki
  void write_asc(ofstream &fp, LINT cid = 0, LINT incrtid=-1);
   static void set_print_cid(BOOLEAN flg=TRUE){print_cid=flg;}

   static void set_tid(LINT ntid){ tid = ntid;}


  LINT size(void) { return nitems; }
};

typedef Transaction *TransactionP;


class CustSeq {
private:
  LINT slen;	// expected number of transactions in sequence
  LINT tlen;	// avg. expected number of items in a transaction
  LINT ntrans;	// number of transactions in sequence
  LINT nitems;	// number of items in sequence
  LINT maxsize;	// size of array trans
  TransactionP *trans;	// transaction in the sequence

public:
  Cid cid;	// customer-id
  CustSeq(Cid cid, LINT seq_len, LINT tot_items);
  ~CustSeq(void);

  BOOLEAN add(String &pat, StringSet &lits);	// adds pattern to transaction
  int write(ofstream &fp);
  int write_asc(ofstream &fp);
  LINT size(void) { return nitems; }
};
