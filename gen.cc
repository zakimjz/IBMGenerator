#include "gen.h"
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <string.h>
#include <math.h>

using namespace std;

//------------------------------- Parameters -------------------------------


void PatternPar::write(ostream &fp)
{
  fp << "\tNumber of patterns = " << npats << endl;
  fp << "\tAverage length of pattern = " << patlen << endl;
  fp << "\tCorrelation between consecutive patterns = " << corr << endl;
  fp << "\tAverage confidence in a rule = " << conf << endl;
  fp << "\tVariation in the confidence = " << conf_var << endl;
}


void TransPar::write(ostream &fp)
{
  fp << "Number of transactions in database = " << ntrans << endl;
  fp << "Average transaction length = " << tlen << endl;
  fp << "Number of items = " << nitems << endl;

  fp << "Large Itemsets:" << endl;
  lits.write(fp);
  fp << endl;
}


// calculate the number of roots, given the number of levels
void TaxPar::calc_values(void)
{
  LINT nset;

  nset = 0;
  nset += (nlevels != 0);
  nset += (fanout != 0);
  nset += (nroots != 0);

  switch (nset)
    {
    case 0:	// fill in defaults
      nroots = 250;
      fanout = 5;
      return;

    case 1:	// need to fill in 1 value
      assert (nlevels == 0);
      if (fanout == 0)
	fanout = 5;
      else if (nroots == 0)
	nroots = 250;
      return;

    case 2:
      if (nlevels == 0)		// all set!
	return;
      if (fanout != 0) {	// calculate nroots
	nroots = (int) (nitems / (1 + pow(DOUBLE(fanout), DOUBLE(nlevels-1))));
	if (nroots < 1)
	  nroots = 1;
      }
      else if (nroots != 0) {	// calculate fanout
	FLOAT temp;
	temp = (FLOAT)nitems / nroots - 1;
	temp = log((DOUBLE)temp) / (nlevels - 1);
	fanout = exp((DOUBLE)temp);
      }
    case 3:			// all set!
      return;
    }
}


void TaxPar::write(ostream &fp)
{
  fp << "Number of transactions in database = " << ntrans << endl;
  fp << "Average transaction length = " << tlen << endl;
  fp << "Number of items = " << nitems << endl;
  fp << "Number of roots = " << nroots << endl;
  fp << "Number of levels = " << nlevels << endl;
  fp << "Average fanout = " << fanout << endl;

  fp << "Large Itemsets:" << endl;
  lits.write(fp);
  fp << endl;
}

//for incremental DB gen
FLOAT SeqPar::CustRatio = 0.0;
FLOAT SeqPar::TransRatio = 0.0;
LINT SeqPar::IncrTid = 0;


void SeqPar::write(ostream &fp)
{
  fp << "Number of customers in database = " << ncust << endl;
  fp << "Average sequence length = " << slen << endl;
  fp << "Average transaction length = " << tlen << endl;
  fp << "Number of items = " << nitems << endl;
  fp << "Repetition-level = " << rept << endl;
  fp << "Variation in repetition-level = " << rept_var << endl;

  fp << "Large Itemsets:" << endl;
  lits.write(fp);
  fp << "Large Sequences:" << endl;
  lseq.write(fp);
  fp << endl;
}


//------------------------------ Taxonomy ------------------------------

const LINT Taxonomy::item_len = 7;

//
// Constructor: reads taxonomy file and builds taxonomy as a DAG
//
Taxonomy::Taxonomy(LINT num_items,	// total number of items
		   LINT num_roots,	// number of roots
		   FLOAT fanout,	// average fanout
		   FLOAT depth_ratio	// average ratio of ....
		   )
  : nitems(num_items), nroots(num_roots), depth(depth_ratio)
{
  LINT i, j;
  LINT next_child;
  PoissonDist nchildren(fanout-1);	// string length

  // allocate memory
  par = new LINT [nitems];
  child_start = new LINT [nitems];
  child_end = new LINT [nitems];

  next_child = nroots;

  // initialize parents (or lack thereof) for roots
  for (i = 0; i < nroots; i++)
    par[i] = -1;

  // set up all the interior nodes
  for (i = 0, j = next_child; i < nitems && next_child < nitems; i++)
    {
      child_start[i] = next_child;
      next_child += nchildren() + 1;
      if (next_child > nitems)
	next_child = nitems;
      child_end[i] = next_child;
      for (; j < next_child; j++)
	par[j] = i;
    }

  // initialize children (or lack thereof) for all the leaves
  for (; i < nitems; i++)
    child_start[i] =
    child_end[i] = -1;
}


Taxonomy::~Taxonomy(void)
{
  delete [] par;
  delete [] child_start;
  delete [] child_end;
}


void Taxonomy::write(ofstream &fp)
{
  for (LINT i = 0; i < nitems; i++)
    if (par[i] >= 0) {
      assert(i != par[i]);
      fp.write((char *)&i, sizeof(LINT));
      fp.write((char *)&par[i], sizeof(LINT));
    }
}


void Taxonomy::write_asc(ofstream &fp)
{
  for (LINT i = 0; i < nitems; i++)
    if (par[i] >= 0) {
      assert(i != par[i]);
      fp << setw(item_len) << i << " " << setw(item_len) << par[i] << endl;
    }
}


void Taxonomy::display(ofstream &fp)
{
  fp << "Taxonomy: " << endl;
  for (LINT i = 0; i < nitems && child_start[i] > 0; i++)
    fp << i << " " << child_start[i] << " " << child_end[i]-1 << endl;
  fp << endl;
}


//------------------------------- ItemSet -------------------------------


ItemSet::ItemSet(LINT num_items, 	// number of items
		 Taxonomy *ptax		// taxonomy (optional)
		 )
  : tax(ptax), nitems(num_items)
{
  ExpDist freq;
  LINT i, j;

  cum_prob = new FLOAT [nitems];
  if (tax)
    tax_prob = new FLOAT [nitems];
  else
    tax_prob = NULL;
  for (i = 0; i < nitems; i++)
    cum_prob[i] = freq();	// prob. that this pattern will be picked

  if (tax) {			// weight(itm) += wieght(children)
    // normalize probabilities for the roots and for children
    normalize(cum_prob, 0, tax->num_roots()-1);
    for (i = 0; i < nitems && tax->num_children(i) > 0; i++)
      normalize(cum_prob, tax->first_child(i), tax->last_child(i));

    // calulate cumulative probabilities for children
    for (i = 0; i < nitems; i++)
      tax_prob[i] = cum_prob[i];
    for (i = 1; i < nitems; i++)
      if (tax->num_children(i) > 0)
	for (j = tax->first_child(i); j < tax->last_child(i); j++)
	  tax_prob[j+1] += tax_prob[j];

    // set real probabilities
    for (i = tax->num_roots(); i < nitems; i++)
      cum_prob[i] *= cum_prob[ tax->parent(i) ] * tax->depth_ratio();
  }

  // normalize probabilites (why -- see get_pat)
  normalize(cum_prob, 0, nitems-1);
  for (i = 1; i < nitems; i++)	// calulate cumulative probabilities
    cum_prob[i] += cum_prob[i-1];
}


ItemSet::~ItemSet()
{
  delete [] cum_prob;
}


//  normalize probabilities between low and high
//
void ItemSet::normalize(FLOAT prob[], LINT low, LINT high)
{
  FLOAT tot;
  LINT i;

  // normalize probabilites
  tot = 0;
  for (i = low; i <= high; i++)
    tot += prob[i];
  for (i = low; i <= high; i++)
    prob[i] /= tot;
}


// returns a pattern chosen at random
//
Item ItemSet::get_item(void)
{
  FLOAT r;
  LINT i;

  // find the desired pattern using cum_prob table
  r = rand();
  // want item i such that cum_prob[i-1] < r <= cum_prob[i];
  i = (int) (r * nitems);			// guess location of item
  i += (int) ((r-cum_prob[i]) * nitems);	// refine guess
  if (i >= nitems)			// check boundaries
    i = nitems-1;
  if (i < 0)
    i = 0;
  while ( i < (nitems-1) && r > cum_prob[i] )	// find item
    i++;
  while ( i > 0 && r <= cum_prob[i-1] )
    i--;
  return i;
}


// if no taxonomy, returns itm
//
Item ItemSet::specialize(Item itm)
{
  FLOAT r;
  LINT i, nchildren;
  Item first, last;

  if (!tax) 		// no taxonomy
    return itm;

  nchildren = tax->num_children(itm);
  if (nchildren == 0)		// no children
    return itm;

  first = tax->child(itm, 0);
  last = tax->child(itm, nchildren-1);

  // find the desired pattern using cum_prob table
  r = rand();
  i = first + (int) (r * nchildren);
  if (i == last)
    i--;
  while ( i < last && r > tax_prob[i] )
    i++;
  while ( i > first && r < tax_prob[i-1] )
    i--;
  return specialize(i);
}


FLOAT ItemSet::weight(Item itm)	// returns prob. of choosing item
{
  if (itm == 0)
    return cum_prob[itm];
  else
    return cum_prob[itm] - cum_prob[itm-1];
}


void ItemSet::display(ofstream &fp)
{
//  if (tax != NULL)
//    tax->display(fp);

  fp << "Items:" << endl;
  fp << setprecision(3);

  if (tax != NULL) {
    if (cum_prob[0] * nitems > 25)
      fp << 0 << "  " << cum_prob[0] * nitems << " "
	<< tax->first_child(0) << " " << tax->last_child(0) << endl;
    for (LINT i = 1; i < nitems; i++)
      if ((cum_prob[i]-cum_prob[i-1]) * nitems > 25)
	fp << i << "  " << (cum_prob[i]-cum_prob[i-1]) * nitems << " "
	  << tax->first_child(i) << " " << tax->last_child(i) << endl;
  }
  else {
    if (cum_prob[0] * nitems > 5)
      fp << 0 << "  " << cum_prob[0] * nitems << endl;
    for (LINT i = 1; i < nitems; i++)
      if ((cum_prob[i]-cum_prob[i-1]) * nitems > 5)
	fp << i << "  " << (cum_prob[i]-cum_prob[i-1]) * nitems << endl;
  }

  fp << setprecision(0);
  fp << endl;
}


//--------------------------------- String ---------------------------------


String::String(LINT n)	// number of items
  : nitems(n)
{
  items = new LINT [nitems];
//  rval = new FLOAT [nitems];
//  ritems = new LINT [nitems];
}


String::~String(void)
{
  delete [] items;
//  delete [] rval;
//  delete [] ritems;
}


void String::display(ofstream &fp, LINT prob_comp)
{
  fp << setw(6) << prob_comp * prob << " " << setw(6) << conf << " ";
  for(LINT i = 0; i < nitems; i++)
    fp << " " << items[i];
  fp << endl;
  return;
}

void String::display(ofstream &fp, StringSet &lits, LINT prob_comp)
{
  LINT i, j;
  StringP lstr;

  fp << setw(6) << prob_comp * prob << " " << setw(6) << conf << " ";
  for(i = 0; i < nitems; i++)
    {
      fp << "  << ";
      lstr = lits.get_pat(items[i]);
      for (j = 0; j < lstr->nitems; j++)
	fp << lstr->items[j] << " ";
      fp << ">>";
    }
  fp << endl;
  return;
}


// shuffles items in string
//
// void String::shuffle(void)
// {
//   UniformDist ud;
//   LINT i, j, ival;
//   FLOAT fval;

//   // associate a random value with each item
//   for (i = 0; i < nitems; i++)
//     rval[i] = ud();

//   // sort items according to the values in rval
//   for (i = 0; i < nitems; i++ )
//     {
//       ival = items[i]; fval = rval[i];
//       for ( j = i; j > 0 && rval[j-1] > fval; j-- ) {
// 	  items[j] = items[j-1];
// 	  rval[j] = rval[j-1];
// 	}
//       items[j] = ival;
//       rval[j] = fval;
//     }
// }


//------------------------------- StringSet -------------------------------


StringSet::StringSet(LINT nitems, 	// number of items
		     PatternPar par,	// npats, patlen, corr, conf & conf_var
		     Taxonomy *ptax,	// taxonomy (optional)
		     FLOAT rept,	// repetition-level
		     FLOAT rept_var	// variation in repetition-level
		     )
  : tax(ptax)
{
  NormalDist conf(par.conf, par.conf_var);
  ExpDist freq;
  ExpDist corr_lvl;
  PoissonDist len(par.patlen-1);	// string length
  NormalDist repeat(rept, rept_var);
  UniformDist ud;

  items = new ItemSet(nitems, tax);	// associate probabilities with items

  LINT i, j, num_same;
  FLOAT tot;

  npats = par.npats;
//  last_pat = 0;
  pat = new StringP [npats];
  for (i = 0; i < npats; i++)
    {
      pat[i] = new String( 1+len() );

      // fill correlated items
      if (par.corr > 0 && i > 0) {	// correlated patterns
	// each pattern has some items same as the previous pattern
	num_same = LINT( pat[i]->size() * par.corr * corr_lvl() + 0.5 );
	if ( num_same > pat[i-1]->size() )
	  num_same = pat[i-1]->size();
	if ( num_same > pat[i]->size() )
	  num_same = pat[i]->size();
	// choose num_same items at random from previous pattern
	Choose shuffle(pat[i-1]->size(), num_same);
	for (j = 0; j < num_same; j++)
	  pat[i]->items[j] = pat[i-1]->item( shuffle.pos(j) );
//	pat[i-1]->shuffle(num_same);
//	for (j = 0; j < num_same; j++)
//	  pat[i]->items[j] = pat[i-1]->rand_item(j);
      }
      else {	// no correlation
	num_same = 0;
      }

      if (rept == 0) { //ERROR: was rept=0 in the oroginal code -- MJZaki
	// fill remaining items at random
	for (j = num_same; j < pat[i]->size(); j++)
	  pat[i]->items[j] = items->get_item();
//	pat[i]->items[j] = LINT(1 + nitems * rand());
      }
      else {
	// some items are repetitions
	FLOAT rept_lvl = repeat();
	for (j = num_same; j < pat[i]->size(); j++)
	  if ( j > 0 && ud() < rept_lvl )	// pick a previous item
	    pat[i]->items[j] = pat[i]->items[ LINT(j*ud()) ];
	  else	// pick random item
	    pat[i]->items[j] = items->get_item();
      }
      pat[i]->prob = freq(); // prob. that this pattern will be picked
      pat[i]->conf = conf(); // used in Transaction::add and CustSeq::add
      			     // to decide how many items to drop from
			     //  this pattern to corrupt it
    }

  if (tax) {
    // weight probabilites with geometric mean of probabilities of items
    for (i = 0; i < npats; i++)
      {
	DOUBLE weight = 1;
	for (j = 0; j < pat[i]->size(); j++)
	  weight *= items->weight(pat[i]->items[j]);
//	cerr << "WEIGHT = " << weight;
	weight = pow(weight, DOUBLE(1)/pat[i]->size());
//	cerr << "  " << weight << endl;
	pat[i]->prob *= weight;
      }
  }

  // normalize probabilites (why -- see get_pat)
  cum_prob = new FLOAT [npats];
  tot = 0;
  for (i = 0; i < npats; i++)
    tot += pat[i]->prob;
  for (i = 0; i < npats; i++)
    pat[i]->prob /= tot;

  // calulate cumulative probabilities
  cum_prob[0] = pat[0]->prob;
  for (i = 1; i < npats; i++)
    cum_prob[i] = cum_prob[i-1] + pat[i]->prob;
//  cerr << cum_prob[npats-1] << endl << flush;

  // allocate space for answer
  LINT maxlen = 0;
  for (i = 1; i < npats; i++)
    if (pat[i]->size() > maxlen)
      maxlen = pat[i]->size();
  answer = new String(maxlen);
}


StringSet::~StringSet()
{
  LINT i;

  for (i = 0; i < npats; i++)
    delete pat[i];
  delete [] pat;
}


// specialize each item in pattern #i and store result in answer
//
StringP StringSet::specialize(LINT i)
{
  answer->set_size( pat[i]->size() );
  answer->set_conf_lvl( pat[i]->conf_lvl() );
  for (LINT j = 0; j < pat[i]->size(); j++)
    answer->set_item(j, items->specialize( pat[i]->item(j) ));
  return answer;
}


// returns pattern #i
//
StringP StringSet::get_pat(LINT i)
{
  if (!tax)
    return pat[i];
  else
    return specialize(i);
}


void StringSet::display(ofstream &fp)
{
  LINT i;

  items->display(fp);

  fp << "ItemSets:" << endl;
  fp << setprecision(3);
  // too lazy to do a sort, so print high-prob. patterns first
  for (i = 0; i < npats; i++)
    if (pat[i]->prob * npats > 10)
      pat[i]->display(fp, npats);
  for (i = 0; i < npats; i++)
    if (pat[i]->prob * npats <= 10 && pat[i]->prob * npats > 1)
      pat[i]->display(fp, npats);
  fp << setprecision(0);
  fp << endl;
}


void StringSet::display(ofstream &fp, StringSet &lits)
{
  LINT i;

  fp << setprecision(3);
  // too lazy to do a sort, so print high-prob. patterns first
  for (i = 0; i < npats; i++)
    if (pat[i]->prob * npats > 6)
      pat[i]->display(fp, lits, npats);
  for (i = 0; i < npats; i++)
    if (pat[i]->prob * npats <= 6)
      pat[i]->display(fp, lits, npats);
  fp << setprecision(0);
}


//------------------------------- StringSet -------------------------------


// returns a pattern chosen at random
//
StringP StringSetIter::get_pat(void)
{
  FLOAT r;
  LINT i = 0;

  if (last_pat < 0) {
    last_pat = -last_pat;
    return strset->pat[last_pat];
  }

  // find the desired pattern using cum_prob table
  r = rand();
  i = (int) (r * strset->npats);
  if (i == strset->npats)
    i--;
  while ( i < (strset->npats-1) && r > strset->cum_prob[i] )
    i++;
  while ( i > 0 && r < strset->cum_prob[i-1] )
    i--;
  last_pat = i;

  if (!strset->tax)
    return strset->pat[i];
  else
    return strset->specialize(i);
}


void StringSetIter::unget_pat(void)
{
  last_pat = -last_pat;
}


//------------------------------ Transaction ------------------------------


// static variables
const LINT Transaction::cid_len = 10;
const LINT Transaction::tid_len = 10;
const LINT Transaction::item_len = 10;
BOOLEAN Transaction::print_cid = TRUE;
LINT Transaction::MAXNITEMS = DEFNITEMS;
LINT Transaction::tid = 0;

Transaction::Transaction(LINT sz)
  : tlen(sz), nitems(0), maxsize(5 * sz)
{
  // maximum size of a transaction is 5 * sz
   if (maxsize > MAXNITEMS) maxsize = MAXNITEMS;
  items = new LINT [maxsize];
}


Transaction::~Transaction()
{
  delete [] items;
}

void Transaction::sort(void)
{
  LINT val;
  LINT i, j;

  for (i = 1; i < nitems; i++ )
    {
      val = items[i];
      for ( j = i; j > 0 && items[j-1] > val; j-- )
	items[j] = items[j-1];
      items[j] = val;
    }
}


BOOLEAN Transaction::add_item(LINT itm)
{
  LINT i;

  for (i = 0; i < nitems; i++)
    if ( items[i] == itm ) return FALSE;

  if (nitems >= maxsize) {	// allocate more memory
    LINT *old_items = items;
    maxsize *= 2;
    items = new LINT [maxsize];
    for (i = 0; i < nitems; i++)
      items[i] = old_items[i];
    delete [] old_items;
  }

  items[nitems++] = itm;
  return TRUE;
}


// adds pattern to transaction
// returns TRUE if pattern was added, FALSE else
//
BOOLEAN Transaction::add(String &pat, BOOLEAN corrupt)
{
  static UniformDist ud;
  LINT i, patlen;

  // corrupt the pattern by reducing its length;
  // conf_lvl for a pattern is decided at the time of pattern creation
  patlen = pat.size();
  if ( corrupt )
    while ( patlen > 0 && ud() > pat.conf_lvl() )
      patlen--;

  // in half of the cases, we drop the pattern that won't fit
  if ( patlen+nitems > tlen )	// not enough space left
    if ( ud() > 0.5 )
      return FALSE;

  // pick "patlen" items at random from pattern
//  if ( patlen < pat.size() )
  Choose shuffle(pat.size(), patlen);
  for (i = 0; i < patlen; i++)
    add_item( pat.item(shuffle.pos(i)) ); // allocates extra space if necessary
//    pat.shuffle(patlen);
//  for (i = 0; i < patlen; i++)
//    add_item( pat.rand_item(i) ); // allocates extra space if necessary

  return TRUE;
}


void Transaction::write(ofstream &fp, LINT cid, LINT incrtid)
{
  if ( nitems == 0 )
    return;
  sort();

  tid++;
  if (cid == 0)		// no customer-id; set cust-id to trans-id
    cid = tid;

  if (print_cid) fp.write((char *)&cid, sizeof(LINT));
  if (incrtid == -1) fp.write((char *)&tid, sizeof(LINT));
  else fp.write((char *)&incrtid, sizeof(LINT));
  fp.write((char *)&nitems, sizeof(LINT));
  fp.write((char *)items, nitems * sizeof(LINT));
  //cout << "TTT " << cid << " " << tid << " " << incrtid << " "
  //<< nitems << endl;
}

void Transaction::write_asc(ofstream &fp, LINT cid, LINT incrtid)
{
  if ( nitems == 0 )
    return;
  sort();

  tid++;
  if (cid == 0)		// no customer-id; set cust-id to trans-id
    cid = tid;

  if (print_cid) fp << cid << " ";
  if (incrtid == -1) fp << tid << " ";
  else  fp << incrtid << " ";
  fp << nitems;
  for (LINT i = 0; i < nitems; i++) {
     fp << " " << items[i];
  }
  fp << endl;
}


//------------------------------ CustSeq ------------------------------


CustSeq::CustSeq(Cid cust_id, LINT seq_len, LINT trans_len)
  : cid(cust_id), slen(seq_len), tlen(trans_len), nitems(0),
    ntrans(seq_len), maxsize(5 * seq_len)
{
  // we reallocate memory if necessary
  trans = new TransactionP [maxsize];
  for (LINT i = 0; i < maxsize; i++)
    trans[i] = NULL;
}


CustSeq::~CustSeq()
{
  for (LINT i = 0; i < maxsize; i++)
    if ( trans[i] )
      delete trans[i];
  delete [] trans;
}


// adds pattern to CustSeq
// returns TRUE if pattern was added, FALSE else
// REWORK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
BOOLEAN CustSeq::add(String &pat, StringSet &lits)
{
  static UniformDist ud;
  LINT i, patlen;
  LINT pos;
  LINT newitems, olditems;
  BOOLEAN corrupt;	// if TRUE, corrupt transactions too

  if ( ud() > pat.conf_lvl() )
    corrupt = TRUE;		// corrupt transactions
  else
    corrupt = FALSE;		// don't corrupt transactions

  // corrupt the pattern by reducing its length;
  // conf_lvl for a pattern is decided at the time of pattern creation
  patlen = pat.size();
  if ( corrupt )
    while ( patlen > 0 && ud() > pat.conf_lvl() )
      patlen--;
  if ( patlen == 0 )	// no trans. left in sequence
    return TRUE;

  // allows transactions to be dropped randomly from the sequence
//  if ( patlen < pat.size() )
  Choose shuffle(pat.size(), patlen);
//    pat.shuffle(patlen);

  // calculate # of items in pattern
  for (newitems = 0, i = 0; i < patlen; i++)
    newitems += lits.get_pat( pat.item( shuffle.pos(i) ) )->size();
//    newitems += lits.get_pat( pat.rand_item(i) )->size();

  // in half of the cases, we drop the pattern that won't fit
  if ( (patlen > slen) || (newitems + nitems > slen * tlen) )
    if ( ud() > 0.5 )
      return FALSE;

  if ( patlen > maxsize ) {	// need to allocate more memory
    TransactionP *old_trans = trans;
    LINT oldsize = maxsize;
    maxsize = patlen*2;
    trans = new TransactionP [maxsize];
    for (i = 0; i < oldsize; i++)
      trans[i] = old_trans[i];
    for (; i < maxsize; i++)
      trans[i] = NULL;
    delete [] old_trans;
  }

  // add new sequence
  Choose *shuffle1 = NULL;
  if (ntrans > patlen)
    shuffle1 = new Choose(ntrans, patlen);
  for (i = 0; i < patlen; i++)
    {
      if ( shuffle1 )
	pos = shuffle1->pos(i);
      else
	pos = i;
      if ( trans[pos] == NULL )
	trans[pos] = new Transaction(tlen);
      olditems = trans[pos]->size();
      trans[pos]->add( *lits.get_pat(pat.item( shuffle.pos(i) )), corrupt );
//      trans[pos]->add( *lits.get_pat(pat.rand_item(i)), corrupt );
      nitems += trans[pos]->size() - olditems;  // update count of #items
    }
  delete shuffle1;

//   pos = ud() * ntrans / patlen;
//   for (i = 0; i < patlen; i++)
//     {
//       if ( trans[pos] == NULL )
// 	trans[pos] = new Transaction(tlen);
//       olditems = trans[pos]->size();
//       trans[pos]->add( *lits.get_pat(pat.item( shuffle.pos(i) )), corrupt );
// //      trans[pos]->add( *lits.get_pat(pat.rand_item(i)), corrupt );
//       nitems += trans[pos]->size() - olditems;  // update count of #items
//       pos += 1 + ud() * ntrans / patlen;
//     }

  return TRUE;
}


int CustSeq::write(ofstream &fp)
{
   LINT cnt = 0; //added by MJZaki

   Transaction::set_tid(0);
   //code to handle incremental DB gen
   LINT incrtid = SeqPar::IncrTid;
   LINT startincrtid = ntrans;
   if (ntrans > 2 && drand48() < SeqPar::CustRatio){
      int ntid = (int) (SeqPar::TransRatio*ntrans);
      if (ntid < 2) ntid = 2;
      startincrtid = ntrans-ntid;
   }

   for (LINT i = 0; i <= ntrans-1; i++)
      if (trans[i] && trans[i]->size() > 0 ){ //added trans[i] -- MJZaki
         if (i >= startincrtid){
            trans[i]->write(fp, cid, incrtid);
            incrtid++;
         }
         else trans[i]->write(fp, cid);
         cnt++; //added by MJZaki
      }
   return cnt; //added by MJZaki
}


int CustSeq::write_asc(ofstream &fp)
{
   LINT cnt=0; //added by MJZaki

   Transaction::set_tid(0);
   //code to handle incremental DB gen
   LINT incrtid = SeqPar::IncrTid;
   LINT startincrtid = ntrans;
   if (ntrans > 2 && drand48() < SeqPar::CustRatio){
      int ntid = (int) (SeqPar::TransRatio*ntrans);
      if (ntid < 2) ntid = 2;
      startincrtid = ntrans-ntid;
   }

   for (LINT i = 0; i <= ntrans-1; i++)
      if (trans[i] &&  trans[i]->size() > 0 ){ //added trans[i] -- MJZaki
         if (i >= startincrtid){
            trans[i]->write_asc(fp, cid, incrtid);
            incrtid++;
         }
         else trans[i]->write_asc(fp, cid);
         cnt++; //added by MJZaki
      }
   return cnt; //added by MJZaki
}
