//
//  VariantContext.cpp
//  OpenGE
//
//  Ported to C++ by Lee Baker on 5/8/12.
//  Open Genomics Engine
//

#include "VariantContext.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
using namespace std;

/**
 * Class VariantContext
 *
 * == High-level overview ==
 *
 * The VariantContext object is a single general class system for representing genetic variation data composed of:
 *
 * * Allele: representing single genetic haplotypes (A, T, ATC, -)
 * * Genotype: an assignment of alleles for each chromosome of a single named sample at a particular locus
 * * VariantContext: an abstract class holding all segregating alleles at a locus as well as genotypes
 *    for multiple individuals containing alleles at that locus
 *
 * The class system works by defining segregating alleles, creating a variant context representing the segregating
 * information at a locus, and potentially creating and associating genotypes with individuals in the context.
 *
 * All of the classes are highly validating -- call validate() if you modify them -- so you can rely on the
 * self-consistency of the data once you have a VariantContext in hand.  The system has a rich set of assessor
 * and manipulator routines, as well as more complex static support routines in VariantContextUtils.
 *
 * The VariantContext (and Genotype) objects are attributed (supporting addition of arbitrary key/value pairs) and
 * filtered (can represent a variation that is viewed as suspect).
 *
 * VariantContexts are dynamically typed, so whether a VariantContext is a SNP, Indel, or NoVariant depends
 * on the properties of the alleles in the context.  See the detailed documentation on the Type parameter below.
 *
 * It's also easy to create subcontexts based on selected genotypes.
 *
 * == Working with Variant Contexts ==
 * By default, VariantContexts are immutable.  In order to access (in the rare circumstances where you need them)
 * setter routines, you need to create MutableVariantContexts and MutableGenotypes.
 *
 * === Some example data ===
 *
 * Allele A, Aref, T, Tref;
 * Allele del, delRef, ATC, ATCref;
 *
 * A [ref] / T at 10
 * GenomeLoc snpLoc = GenomeLocParser.createGenomeLoc("chr1", 10, 10);
 *
 * - / ATC [ref] from 20-23
 * GenomeLoc delLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 22);
 *
 *  // - [ref] / ATC immediately after 20
 * GenomeLoc insLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 20);
 *
 * === Alleles ===
 *
 * See the documentation in the Allele class itself
 *
 * What are they?
 *
 * Alleles can be either reference or non-reference
 *
 * Example alleles used here:
 *
 *   del = new Allele("-");
 *   A = new Allele("A");
 *   Aref = new Allele("A", true);
 *   T = new Allele("T");
 *   ATC = new Allele("ATC");
 *
 * === Creating variant contexts ===
 *
 * ==== By hand ====
 *
 * Here's an example of a A/T polymorphism with the A being reference:
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, snpLoc, Arrays.asList(Aref, T));
 * </pre>
 *
 * If you want to create a non-variant site, just put in a single reference allele
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, snpLoc, Arrays.asList(Aref));
 * </pre>
 *
 * A deletion is just as easy:
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, delLoc, Arrays.asList(ATCref, del));
 * </pre>
 *
 * The only 2 things that distinguishes between a insertion and deletion are the reference allele
 * and the location of the variation.  An insertion has a Null reference allele and at least
 * one non-reference Non-Null allele.  Additionally, the location of the insertion is immediately after
 * a 1-bp GenomeLoc (at say 20).
 *
 * <pre>
 * VariantContext vc = new VariantContext("name", insLoc, Arrays.asList(delRef, ATC));
 * </pre>
 *
 * ==== Converting rods and other data structures to VCs ====
 *
 * You can convert many common types into VariantContexts using the general function:
 *
 * <pre>
 * VariantContextAdaptors.convertToVariantContext(name, myObject)
 * </pre>
 *
 * dbSNP and VCFs, for example, can be passed in as myObject and a VariantContext corresponding to that
 * object will be returned.  A null return type indicates that the type isn't yet supported.  This is the best
 * and easiest way to create contexts using RODs.
 *
 *
 * === Working with genotypes ===
 *
 * <pre>
 * List<Allele> alleles = Arrays.asList(Aref, T);
 * Genotype g1 = new Genotype(Arrays.asList(Aref, Aref), "g1", 10);
 * Genotype g2 = new Genotype(Arrays.asList(Aref, T), "g2", 10);
 * Genotype g3 = new Genotype(Arrays.asList(T, T), "g3", 10);
 * VariantContext vc = new VariantContext(snpLoc, alleles, Arrays.asList(g1, g2, g3));
 * </pre>
 *
 * At this point we have 3 genotypes in our context, g1-g3.
 *
 * You can assess a good deal of information about the genotypes through the VariantContext:
 *
 * <pre>
 * vc.hasGenotypes()
 * vc.isMonomorphicInSamples()
 * vc.isPolymorphicInSamples()
 * vc.getSamples().size()
 *
 * vc.getGenotypes()
 * vc.getGenotypes().get("g1")
 * vc.hasGenotype("g1")
 *
 * vc.getCalledChrCount()
 * vc.getCalledChrCount(Aref)
 * vc.getCalledChrCount(T)
 * </pre>
 *
 * === NO_CALL alleles ===
 *
 * The system allows one to create Genotypes carrying special NO_CALL alleles that aren't present in the
 * set of context alleles and that represent undetermined alleles in a genotype:
 *
 * Genotype g4 = new Genotype(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), "NO_DATA_FOR_SAMPLE", 10);
 *
 *
 * === subcontexts ===
 * It's also very easy get subcontext based only the data in a subset of the genotypes:
 *
 * <pre>
 * VariantContext vc12 = vc.subContextFromGenotypes(Arrays.asList(g1,g2));
 * VariantContext vc1 = vc.subContextFromGenotypes(Arrays.asList(g1));
 * </pre>
 *
 * @author depristo
 */

#if 0
class VariantContext {
protected CommonInfo commonInfo = null;
public final static double NO_LOG10_PERROR = CommonInfo.NO_LOG10_PERROR;

@Deprecated // ID is no longer stored in the attributes map
private final static String ID_KEY = "ID";

private final Byte REFERENCE_BASE_FOR_INDEL;

public final static Set<String> PASSES_FILTERS = Collections.unmodifiableSet(new LinkedHashSet<String>());

/** The location of this VariantContext */
private final String ID;

/** A mapping from sampleName -> genotype objects for all genotypes associated with this context */
protected GenotypesContext genotypes = null;

/** Counts for each of the possible Genotype types in this context */
protected int[] genotypeCounts = null;

public final static GenotypesContext NO_GENOTYPES = GenotypesContext.NO_GENOTYPES;



/* cached monomorphic value: null -> not yet computed, False, True */
private Boolean monomorphic = null;

// ---------------------------------------------------------------------------------------------------------
//
// validation mode
//
// ---------------------------------------------------------------------------------------------------------

public enum Validation {
    REF_PADDING,
    ALLELES,
    GENOTYPES
}

private final static EnumSet<Validation> ALL_VALIDATION = EnumSet.allOf(Validation.class);
private final static EnumSet<Validation> NO_VALIDATION = EnumSet.noneOf(Validation.class);

// ---------------------------------------------------------------------------------------------------------
//
// constructors: see VariantContextBuilder
//
// ---------------------------------------------------------------------------------------------------------

/**
 * Copy constructor
 *
 * @param other the VariantContext to copy
 */
protected VariantContext(VariantContext other) {
    this(other.getSource(), other.getID(), other.getChr(), other.getStart(), other.getEnd(),
         other.getAlleles(), other.getGenotypes(), other.getLog10PError(),
         other.getFiltersMaybeNull(),
         other.getAttributes(), other.REFERENCE_BASE_FOR_INDEL,
         NO_VALIDATION);
}

/**
 * the actual constructor.  Private access only
 *
 * @param source          source
 * @param contig          the contig
 * @param start           the start base (one based)
 * @param stop            the stop reference base (one based)
 * @param alleles         alleles
 * @param genotypes       genotypes map
 * @param log10PError  qual
 * @param filters         filters: use null for unfiltered and empty set for passes filters
 * @param attributes      attributes
 * @param referenceBaseForIndel   padded reference base
 * @param validationToPerform     set of validation steps to take
 */
protected VariantContext(String source, String ID,
                         String contig, long start, long stop,
                         Collection<Allele> alleles, GenotypesContext genotypes,
                         double log10PError, Set<String> filters, Map<String, Object> attributes,
                         Byte referenceBaseForIndel,
                         EnumSet<Validation> validationToPerform ) {
    if ( contig == null ) { throw new IllegalArgumentException("Contig cannot be null"); }
    this.contig = contig;
    this.start = start;
    this.stop = stop;
    
    // intern for efficiency.  equals calls will generate NPE if ID is inappropriately passed in as null
    if ( ID == null || ID.equals("") ) throw new IllegalArgumentException("ID field cannot be the null or the empty string");
    this.ID = ID.equals(VCFConstants.EMPTY_ID_FIELD) ? VCFConstants.EMPTY_ID_FIELD : ID;
    
    this.commonInfo = new CommonInfo(source, log10PError, filters, attributes);
    REFERENCE_BASE_FOR_INDEL = referenceBaseForIndel;
    
    // todo -- remove me when this check is no longer necessary
    if ( this.commonInfo.hasAttribute(ID_KEY) )
        throw new IllegalArgumentException("Trying to create a VariantContext with a ID key.  Please use provided constructor argument ID");
    
    if ( alleles == null ) { throw new IllegalArgumentException("Alleles cannot be null"); }
    
    // we need to make this a LinkedHashSet in case the user prefers a given ordering of alleles
    this.alleles = makeAlleles(alleles);
    
    if ( genotypes == null || genotypes == NO_GENOTYPES ) {
        this.genotypes = NO_GENOTYPES;
    } else {
        this.genotypes = genotypes.immutable();
    }
    
    // cache the REF and ALT alleles
    int nAlleles = alleles.size();
    for ( Allele a : alleles ) {
        if ( a.isReference() ) {
            REF = a;
        } else if ( nAlleles == 2 ) { // only cache ALT when biallelic
            ALT = a;
        }
    }
    
    if ( ! validationToPerform.isEmpty() ) {
        validate(validationToPerform);
    }
}

// ---------------------------------------------------------------------------------------------------------
//
// Selectors
//
// ---------------------------------------------------------------------------------------------------------

public VariantContext subContextFromSamples(Set<String> sampleNames, Collection<Allele> alleles) {
    VariantContextBuilder builder = new VariantContextBuilder(this);
    return builder.genotypes(genotypes.subsetToSamples(sampleNames)).alleles(alleles).make();
}

public VariantContext subContextFromSamples(Set<String> sampleNames) {
    VariantContextBuilder builder = new VariantContextBuilder(this);
    GenotypesContext newGenotypes = genotypes.subsetToSamples(sampleNames);
    return builder.genotypes(newGenotypes).alleles(allelesOfGenotypes(newGenotypes)).make();
}

public VariantContext subContextFromSample(String sampleName) {
    return subContextFromSamples(Collections.singleton(sampleName));
}

/**
 * helper routine for subcontext
 * @param genotypes genotypes
 * @return allele set
 */
private final Set<Allele> allelesOfGenotypes(Collection<Genotype> genotypes) {
    final Set<Allele> alleles = new HashSet<Allele>();
    
    boolean addedref = false;
    for ( final Genotype g : genotypes ) {
        for ( final Allele a : g.getAlleles() ) {
            addedref = addedref || a.isReference();
            if ( a.isCalled() )
                alleles.add(a);
        }
    }
    if ( ! addedref ) alleles.add(getReference());
    
    return alleles;
}

// ---------------------------------------------------------------------------------------------------------
//
// type operations
//
// ---------------------------------------------------------------------------------------------------------

/**
 * see: http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=handbook&part=ch5&rendertype=table&id=ch5.ch5_t3
 *
 * Format:
 * dbSNP variation class
 * Rules for assigning allele classes
 * Sample allele definition
 *
 * Single Nucleotide Polymorphisms (SNPs)a
 *   Strictly defined as single base substitutions involving A, T, C, or G.
 *   A/T
 *
 * Deletion/Insertion Polymorphisms (DIPs)
 *   Designated using the full sequence of the insertion as one allele, and either a fully
 *   defined string for the variant allele or a '-' character to specify the deleted allele.
 *   This class will be assigned to a variation if the variation alleles are of different lengths or
 *   if one of the alleles is deleted ('-').
 *   T/-/CCTA/G
 *
 * No-variation
 *   Reports may be submitted for segments of sequence that are assayed and determined to be invariant
 *   in the sample.
 *   (NoVariation)
 *
 * Mixed
 *   Mix of other classes
 *
 * Also supports NO_VARIATION type, used to indicate that the site isn't polymorphic in the population
 *
 *
 * Not currently supported:
 *
 * Heterozygous sequencea
 * The term heterozygous is used to specify a region detected by certain methods that do not
 * resolve the polymorphism into a specific sequence motif. In these cases, a unique flanking
 * sequence must be provided to define a sequence context for the variation.
 * (heterozygous)
 *
 * Microsatellite or short tandem repeat (STR)
 * Alleles are designated by providing the repeat motif and the copy number for each allele.
 * Expansion of the allele repeat motif designated in dbSNP into full-length sequence will
 * be only an approximation of the true genomic sequence because many microsatellite markers are
 * not fully sequenced and are resolved as size variants only.
 * (CAC)8/9/10/11
 *
 * Named variant
 * Applies to insertion/deletion polymorphisms of longer sequence features, such as retroposon
 * dimorphism for Alu or line elements. These variations frequently include a deletion '-' indicator
 * for the absent allele.
 * (alu) / -
 *
 * Multi-Nucleotide Polymorphism (MNP)
 *   Assigned to variations that are multi-base variations of a single, common length
 *   GGA/AGT
 */
public enum Type {
    NO_VARIATION,
    SNP,
    MNP,    // a multi-nucleotide polymorphism
    INDEL,
    SYMBOLIC,
    MIXED,
}
#endif
/**
 * Determines (if necessary) and returns the type of this variation by examining the alleles it contains.
 *
 * @return the type of this VariantContext
 **/
VariantContext::Type VariantContext::getType() {
    if ( type == NULL )
        determineType();
    
    return *type;
}
#if 0

/**
 * convenience method for SNPs
 *
 * @return true if this is a SNP, false otherwise
 */
public boolean isSNP() { return getType() == Type.SNP; }


/**
 * convenience method for variants
 *
 * @return true if this is a variant allele, false if it's reference
 */
public boolean isVariant() { return getType() != Type.NO_VARIATION; }

/**
 * convenience method for point events
 *
 * @return true if this is a SNP or ref site, false if it's an indel or mixed event
 */
public boolean isPointEvent() { return isSNP() || !isVariant(); }
    
#endif
/**
 * convenience method for indels
 *
 * @return true if this is an indel, false otherwise
 */
bool VariantContext::isIndel() { return getType() == INDEL; }


/**
 * @return true if the alleles indicate a simple insertion (i.e., the reference allele is Null)
 */
bool VariantContext::isSimpleInsertion() {
    // can't just call !isSimpleDeletion() because of complex indels
    return getType() == VariantContext::INDEL && getReference().isNull() && isBiallelic();
}

/**
 * @return true if the alleles indicate a simple deletion (i.e., a single alt allele that is Null)
 */
bool VariantContext::isSimpleDeletion() {
    // can't just call !isSimpleInsertion() because of complex indels
    return getType() == INDEL && getAlternateAllele(0).isNull() && isBiallelic();
}

/**
 * @return true if the alleles indicate neither a simple deletion nor a simple insertion
 */
bool VariantContext::isComplexIndel() {
    return isIndel() && !isSimpleDeletion() && !isSimpleInsertion();
}

bool VariantContext::isSymbolic() {
    return getType() == SYMBOLIC;
}

bool VariantContext::isMNP() {
    return getType() == MNP;
}
#if 0

/**
 * convenience method for indels
 *
 * @return true if this is an mixed variation, false otherwise
 */
public boolean isMixed() { return getType() == Type.MIXED; }


// ---------------------------------------------------------------------------------------------------------
//
// Generic accessors
//
// ---------------------------------------------------------------------------------------------------------

public boolean hasID() {
    return getID() != VCFConstants.EMPTY_ID_FIELD;
}

public boolean emptyID() {
    return ! hasID();
}

public String getID() {
    return ID;
}

public boolean hasReferenceBaseForIndel() {
    return REFERENCE_BASE_FOR_INDEL != null;
}

// the indel base that gets stripped off for indels
public Byte getReferenceBaseForIndel() {
    return REFERENCE_BASE_FOR_INDEL;
}

// ---------------------------------------------------------------------------------------------------------
//
// get routines to access context info fields
//
// ---------------------------------------------------------------------------------------------------------
public String getSource()                   { return commonInfo.getName(); }
public Set<String> getFiltersMaybeNull()    { return commonInfo.getFiltersMaybeNull(); }
public Set<String> getFilters()             { return commonInfo.getFilters(); }
public boolean isFiltered()                 { return commonInfo.isFiltered(); }
public boolean isNotFiltered()              { return commonInfo.isNotFiltered(); }
public boolean filtersWereApplied()         { return commonInfo.filtersWereApplied(); }
public boolean hasLog10PError()             { return commonInfo.hasLog10PError(); }
public double getLog10PError()              { return commonInfo.getLog10PError(); }
public double getPhredScaledQual()          { return commonInfo.getPhredScaledQual(); }

public Map<String, Object>  getAttributes() { return commonInfo.getAttributes(); }
public boolean hasAttribute(String key)     { return commonInfo.hasAttribute(key); }
public Object getAttribute(String key)      { return commonInfo.getAttribute(key); }

public Object getAttribute(String key, Object defaultValue) {
    return commonInfo.getAttribute(key, defaultValue);
}

public String getAttributeAsString(String key, String defaultValue)   { return commonInfo.getAttributeAsString(key, defaultValue); }
public int getAttributeAsInt(String key, int defaultValue)            { return commonInfo.getAttributeAsInt(key, defaultValue); }
public double getAttributeAsDouble(String key, double  defaultValue)  { return commonInfo.getAttributeAsDouble(key, defaultValue); }
public boolean getAttributeAsBoolean(String key, boolean  defaultValue)  { return commonInfo.getAttributeAsBoolean(key, defaultValue); }

// ---------------------------------------------------------------------------------------------------------
//
// Working with alleles
//
// ---------------------------------------------------------------------------------------------------------
#endif
/**
 * @return the reference allele for this context
 */
Allele VariantContext::getReference() {
    Allele * ref = REF;
    assert(ref != NULL);
    if ( ref == NULL ) {
        cerr << "BUG: no reference allele found at " << endl;
        abort();
    }
    return *ref;
}


/**
 * @return true if the context is strictly bi-allelic
 */
bool VariantContext::isBiallelic() {
    return getNAlleles() == 2;
}

/**
 * @return The number of segregating alleles in this context
 */
int VariantContext::getNAlleles() {
    return alleles.size();
}
#if 0

/**
 * @return The allele sharing the same bases as this String.  A convenience method; better to use byte[]
 */
public Allele getAllele(String allele) {
    return getAllele(allele.getBytes());
}

/**
 * @return The allele sharing the same bases as this byte[], or null if no such allele is present.
 */
public Allele getAllele(byte[] allele) {
    return Allele.getMatchingAllele(getAlleles(), allele);
}

/**
 * @return True if this context contains Allele allele, or false otherwise
 */
public boolean hasAllele(Allele allele) {
    return hasAllele(allele, false);
}

public boolean hasAllele(Allele allele, boolean ignoreRefState) {
    if ( allele == REF || allele == ALT ) // optimization for cached cases
        return true;
    
    for ( Allele a : getAlleles() ) {
        if ( a.equals(allele, ignoreRefState) )
            return true;
    }
    
    return false;
}


/**
 * Gets the alleles.  This method should return all of the alleles present at the location,
 * including the reference allele.  There are no constraints imposed on the ordering of alleles
 * in the set. If the reference is not an allele in this context it will not be included.
 *
 * @return the set of alleles
 */
public List<Allele> getAlleles() { return alleles; }
#endif
/**
 * Gets the alternate alleles.  This method should return all the alleles present at the location,
 * NOT including the reference allele.  There are no constraints imposed on the ordering of alleles
 * in the set.
 *
 * @return the set of alternate alleles
 */
vector<Allele> VariantContext::getAlternateAlleles() {
    return vector<Allele>(alleles.begin() + 1, alleles.end());
}
#if 0

/**
 * Gets the sizes of the alternate alleles if they are insertion/deletion events, and returns a list of their sizes
 *
 * @return a list of indel lengths ( null if not of type indel or mixed )
 */
public List<Integer> getIndelLengths() {
    if ( getType() != Type.INDEL && getType() != Type.MIXED ) {
        return null;
    }
    
    List<Integer> lengths = new ArrayList<Integer>();
    for ( Allele a : getAlternateAlleles() ) {
        lengths.add(a.length() - getReference().length());
    }
    
    return lengths;
}
#endif

/**
 * @param i -- the ith allele (from 0 to n - 2 for a context with n alleles including a reference allele)
 * @return the ith non-reference allele in this context
 * @throws IllegalArgumentException if i is invalid
 */
Allele VariantContext::getAlternateAllele(int i) {
    return alleles[i+1];
}
#if 0

/**
 * @param  other  VariantContext whose alternate alleles to compare against
 * @return true if this VariantContext has the same alternate alleles as other,
 *         regardless of ordering. Otherwise returns false.
 */
public boolean hasSameAlternateAllelesAs ( VariantContext other ) {
    List<Allele> thisAlternateAlleles = getAlternateAlleles();
    List<Allele> otherAlternateAlleles = other.getAlternateAlleles();
    
    if ( thisAlternateAlleles.size() != otherAlternateAlleles.size() ) {
        return false;
    }
    
    for ( Allele allele : thisAlternateAlleles ) {
        if ( ! otherAlternateAlleles.contains(allele) ) {
            return false;
        }
    }
    
    return true;
}

// ---------------------------------------------------------------------------------------------------------
//
// Working with genotypes
//
// ---------------------------------------------------------------------------------------------------------

/**
 * @return the number of samples in the context
 */
public int getNSamples() {
    return genotypes.size();
}

/**
 * @return true if the context has associated genotypes
 */
public boolean hasGenotypes() {
    return ! genotypes.isEmpty();
}

public boolean hasGenotypes(Collection<String> sampleNames) {
    return genotypes.containsSamples(sampleNames);
}

/**
 * @return set of all Genotypes associated with this context
 */
public GenotypesContext getGenotypes() {
    return genotypes;
}

public Iterable<Genotype> getGenotypesOrderedByName() {
    return genotypes.iterateInSampleNameOrder();
}

public Iterable<Genotype> getGenotypesOrderedBy(Iterable<String> sampleOrdering) {
    return genotypes.iterateInSampleNameOrder(sampleOrdering);
}

/**
 * Returns a map from sampleName -> Genotype for the genotype associated with sampleName.  Returns a map
 * for consistency with the multi-get function.
 *
 * @param sampleName
 * @return
 * @throws IllegalArgumentException if sampleName isn't bound to a genotype
 */
public GenotypesContext getGenotypes(String sampleName) {
    return getGenotypes(Collections.singleton(sampleName));
}

/**
 * Returns a map from sampleName -> Genotype for each sampleName in sampleNames.  Returns a map
 * for consistency with the multi-get function.
 *
 * For testing convenience only
 *
 * @param sampleNames a unique list of sample names
 * @return
 * @throws IllegalArgumentException if sampleName isn't bound to a genotype
 */
protected GenotypesContext getGenotypes(Collection<String> sampleNames) {
    return getGenotypes().subsetToSamples(new HashSet<String>(sampleNames));
}

public GenotypesContext getGenotypes(Set<String> sampleNames) {
    return getGenotypes().subsetToSamples(sampleNames);
}


/**
 * @return the set of all sample names in this context, not ordered
 */
public Set<String> getSampleNames() {
    return getGenotypes().getSampleNames();
}

public List<String> getSampleNamesOrderedByName() {
    return getGenotypes().getSampleNamesOrderedByName();
}

/**
 * @param sample  the sample name
 *
 * @return the Genotype associated with the given sample in this context or null if the sample is not in this context
 */
public Genotype getGenotype(String sample) {
    return getGenotypes().get(sample);
}

public boolean hasGenotype(String sample) {
    return getGenotypes().containsSample(sample);
}

public Genotype getGenotype(int ith) {
    return genotypes.get(ith);
}


/**
 * Returns the number of chromosomes carrying any allele in the genotypes (i.e., excluding NO_CALLS)
 *
 * @return chromosome count
 */
public int getCalledChrCount() {
    int n = 0;
    
    for ( final Genotype g : getGenotypes() ) {
        for ( final Allele a : g.getAlleles() )
            n += a.isNoCall() ? 0 : 1;
    }
    
    return n;
}

/**
 * Returns the number of chromosomes carrying allele A in the genotypes
 *
 * @param a allele
 * @return chromosome count
 */
public int getCalledChrCount(Allele a) {
    int n = 0;
    
    for ( final Genotype g : getGenotypes() ) {
        n += g.getAlleles(a).size();
    }
    
    return n;
}

/**
 * Genotype-specific functions -- are the genotypes monomorphic w.r.t. to the alleles segregating at this
 * site?  That is, is the number of alternate alleles among all fo the genotype == 0?
 *
 * @return true if it's monomorphic
 */
public boolean isMonomorphicInSamples() {
    if ( monomorphic == null )
        monomorphic = ! isVariant() || (hasGenotypes() && getCalledChrCount(getReference()) == getCalledChrCount());
    return monomorphic;
}

/**
 * Genotype-specific functions -- are the genotypes polymorphic w.r.t. to the alleles segregating at this
 * site?  That is, is the number of alternate alleles among all fo the genotype > 0?
 *
 * @return true if it's polymorphic
 */
public boolean isPolymorphicInSamples() {
    return ! isMonomorphicInSamples();
}

private void calculateGenotypeCounts() {
    if ( genotypeCounts == null ) {
        genotypeCounts = new int[Genotype.Type.values().length];
        
        for ( final Genotype g : getGenotypes() ) {
            genotypeCounts[g.getType().ordinal()]++;
        }
    }
}

/**
 * Genotype-specific functions -- how many no-calls are there in the genotypes?
 *
 * @return number of no calls
 */
public int getNoCallCount() {
    calculateGenotypeCounts();
    return genotypeCounts[Genotype.Type.NO_CALL.ordinal()];
}

/**
 * Genotype-specific functions -- how many hom ref calls are there in the genotypes?
 *
 * @return number of hom ref calls
 */
public int getHomRefCount() {
    calculateGenotypeCounts();
    return genotypeCounts[Genotype.Type.HOM_REF.ordinal()];
}

/**
 * Genotype-specific functions -- how many het calls are there in the genotypes?
 *
 * @return number of het calls
 */
public int getHetCount() {
    calculateGenotypeCounts();
    return genotypeCounts[Genotype.Type.HET.ordinal()];
}

/**
 * Genotype-specific functions -- how many hom var calls are there in the genotypes?
 *
 * @return number of hom var calls
 */
public int getHomVarCount() {
    return genotypeCounts[Genotype.Type.HOM_VAR.ordinal()];
}

/**
 * Genotype-specific functions -- how many mixed calls are there in the genotypes?
 *
 * @return number of mixed calls
 */
public int getMixedCount() {
    return genotypeCounts[Genotype.Type.MIXED.ordinal()];
}

// ---------------------------------------------------------------------------------------------------------
//
// validation: extra-strict validation routines for paranoid users
//
// ---------------------------------------------------------------------------------------------------------

/**
 * Run all extra-strict validation tests on a Variant Context object
 *
 * @param reference        the true reference allele
 * @param paddedRefBase    the reference base used for padding indels
 * @param rsIDs            the true dbSNP IDs
 */
public void extraStrictValidation(Allele reference, Byte paddedRefBase, Set<String> rsIDs) {
    // validate the reference
    validateReferenceBases(reference, paddedRefBase);
    
    // validate the RS IDs
    validateRSIDs(rsIDs);
    
    // validate the altenate alleles
    validateAlternateAlleles();
    
    // validate the AN and AC fields
    validateChromosomeCounts();
    
    // TODO: implement me
    //checkReferenceTrack();
}

public void validateReferenceBases(Allele reference, Byte paddedRefBase) {
    if ( reference == null )
        return;
    
    // don't validate if we're a complex event
    if ( !isComplexIndel() && !reference.isNull() && !reference.basesMatch(getReference()) ) {
        throw new TribbleException.InternalCodecException(String.format("the REF allele is incorrect for the record at position %s:%d, fasta says %s vs. VCF says %s", getChr(), getStart(), reference.getBaseString(), getReference().getBaseString()));
    }
    
    // we also need to validate the padding base for simple indels
    if ( hasReferenceBaseForIndel() && !getReferenceBaseForIndel().equals(paddedRefBase) ) {
        throw new TribbleException.InternalCodecException(String.format("the padded REF base is incorrect for the record at position %s:%d, fasta says %s vs. VCF says %s", getChr(), getStart(), (char)paddedRefBase.byteValue(), (char)getReferenceBaseForIndel().byteValue()));
    }
}

public void validateRSIDs(Set<String> rsIDs) {
    if ( rsIDs != null && hasID() ) {
        for ( String id : getID().split(VCFConstants.ID_FIELD_SEPARATOR) ) {
            if ( id.startsWith("rs") && !rsIDs.contains(id) )
                throw new TribbleException.InternalCodecException(String.format("the rsID %s for the record at position %s:%d is not in dbSNP", id, getChr(), getStart()));
        }
    }
}

public void validateAlternateAlleles() {
    if ( !hasGenotypes() )
        return;
    
    List<Allele> reportedAlleles = getAlleles();
    Set<Allele> observedAlleles = new HashSet<Allele>();
    observedAlleles.add(getReference());
    for ( final Genotype g : getGenotypes() ) {
        if ( g.isCalled() )
            observedAlleles.addAll(g.getAlleles());
    }
    
    if ( reportedAlleles.size() != observedAlleles.size() )
        throw new TribbleException.InternalCodecException(String.format("the ALT allele(s) for the record at position %s:%d do not match what is observed in the per-sample genotypes", getChr(), getStart()));
    
    int originalSize = reportedAlleles.size();
    // take the intersection and see if things change
    observedAlleles.retainAll(reportedAlleles);
    if ( observedAlleles.size() != originalSize )
        throw new TribbleException.InternalCodecException(String.format("the ALT allele(s) for the record at position %s:%d do not match what is observed in the per-sample genotypes", getChr(), getStart()));
}

public void validateChromosomeCounts() {
    if ( !hasGenotypes() )
        return;
    
    // AN
    if ( hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) ) {
        int reportedAN = Integer.valueOf(getAttribute(VCFConstants.ALLELE_NUMBER_KEY).toString());
        int observedAN = getCalledChrCount();
        if ( reportedAN != observedAN )
            throw new TribbleException.InternalCodecException(String.format("the Allele Number (AN) tag is incorrect for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedAN, observedAN));
    }
    
    // AC
    if ( hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
        ArrayList<Integer> observedACs = new ArrayList<Integer>();
        
        // if there are alternate alleles, record the relevant tags
        if ( getAlternateAlleles().size() > 0 ) {
            for ( Allele allele : getAlternateAlleles() ) {
                observedACs.add(getCalledChrCount(allele));
            }
        }
        else { // otherwise, set them to 0
            observedACs.add(0);
        }
        
        if ( getAttribute(VCFConstants.ALLELE_COUNT_KEY) instanceof List ) {
            Collections.sort(observedACs);
            List reportedACs = (List)getAttribute(VCFConstants.ALLELE_COUNT_KEY);
            Collections.sort(reportedACs);
            if ( observedACs.size() != reportedACs.size() )
                throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag doesn't have the correct number of values for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedACs.size(), observedACs.size()));
            for (int i = 0; i < observedACs.size(); i++) {
                if ( Integer.valueOf(reportedACs.get(i).toString()) != observedACs.get(i) )
                    throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag is incorrect for the record at position %s:%d, %s vs. %d", getChr(), getStart(), reportedACs.get(i), observedACs.get(i)));
            }
        } else {
            if ( observedACs.size() != 1 )
                throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag doesn't have enough values for the record at position %s:%d", getChr(), getStart()));
            int reportedAC = Integer.valueOf(getAttribute(VCFConstants.ALLELE_COUNT_KEY).toString());
            if ( reportedAC != observedACs.get(0) )
                throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag is incorrect for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedAC, observedACs.get(0)));
        }
    }
}

// ---------------------------------------------------------------------------------------------------------
//
// validation: the normal validation routines are called automatically upon creation of the VC
//
// ---------------------------------------------------------------------------------------------------------

private boolean validate(final EnumSet<Validation> validationToPerform) {
    for (final Validation val : validationToPerform ) {
        switch (val) {
            case ALLELES: validateAlleles(); break;
            case REF_PADDING: validateReferencePadding(); break;
            case GENOTYPES: validateGenotypes(); break;
            default: throw new IllegalArgumentException("Unexpected validation mode " + val);
        }
    }
    
    return true;
}

private void validateReferencePadding() {
    if (hasSymbolicAlleles()) // symbolic alleles don't need padding...
        return;
    
    boolean needsPadding = (getReference().length() == getEnd() - getStart()); // off by one because padded base was removed
    
    if ( needsPadding && !hasReferenceBaseForIndel() )
        throw new ReviewedStingException("Badly formed variant context at location " + getChr() + ":" + getStart() + "; no padded reference base was provided.");
}

private void validateAlleles() {
    // check alleles
    boolean alreadySeenRef = false, alreadySeenNull = false;
    for ( Allele allele : alleles ) {
        // make sure there's only one reference allele
        if ( allele.isReference() ) {
            if ( alreadySeenRef ) throw new IllegalArgumentException("BUG: Received two reference tagged alleles in VariantContext " + alleles + " this=" + this);
            alreadySeenRef = true;
        }
        
        if ( allele.isNoCall() ) {
            throw new IllegalArgumentException("BUG: Cannot add a no call allele to a variant context " + alleles + " this=" + this);
        }
        
        // make sure there's only one null allele
        if ( allele.isNull() ) {
            if ( alreadySeenNull ) throw new IllegalArgumentException("BUG: Received two null alleles in VariantContext " + alleles + " this=" + this);
            alreadySeenNull = true;
        }
    }
    
    // make sure there's one reference allele
    if ( ! alreadySeenRef )
        throw new IllegalArgumentException("No reference allele found in VariantContext");
    
    //        if ( getType() == Type.INDEL ) {
    //            if ( getReference().length() != (getLocation().size()-1) ) {
    long length = (stop - start) + 1;
    if ( (getReference().isNull() && length != 1 ) ||
        (getReference().isNonNull() && (length - getReference().length()  > 1))) {
        throw new IllegalStateException("BUG: GenomeLoc " + contig + ":" + start + "-" + stop + " has a size == " + length + " but the variation reference allele has length " + getReference().length() + " this = " + this);
    }
}

private void validateGenotypes() {
    if ( this.genotypes == null ) throw new IllegalStateException("Genotypes is null");
    
    for ( final Genotype g : this.genotypes ) {
        if ( g.isAvailable() ) {
            for ( Allele gAllele : g.getAlleles() ) {
                if ( ! hasAllele(gAllele) && gAllele.isCalled() )
                    throw new IllegalStateException("Allele in genotype " + gAllele + " not in the variant context " + alleles);
            }
        }
    }
}
#endif

// ---------------------------------------------------------------------------------------------------------
//
// utility routines
//
// ---------------------------------------------------------------------------------------------------------

void VariantContext::determineType() {
    if ( type == NULL ) {
        switch ( getNAlleles() ) {
            case 0:
                cerr << "Unexpected error: requested type of VariantContext with no alleles!" << endl;
                abort();
            case 1:
                // note that this doesn't require a reference allele.  You can be monomorphic independent of having a
                // reference allele
                type = new Type(NO_VARIATION);
                break;
            default:
                determinePolymorphicType();
        }
    }
}

void VariantContext::determinePolymorphicType() {
    type = NULL;
    
    // do a pairwise comparison of all alleles against the reference allele
    for ( vector<Allele>::iterator allele = alleles.begin() ; allele != alleles.end(); allele++) {
        if (*allele == *REF )
            continue;
        
        // find the type of this allele relative to the reference
        Type biallelicType = typeOfBiallelicVariant(*REF, *allele);
        
        // for the first alternate allele, set the type to be that one
        if ( type == NULL ) {
            type = new Type(biallelicType);
        }
        // if the type of this allele is different from that of a previous one, assign it the MIXED type and quit
        else if ( biallelicType != *type ) {
            type = new Type(MIXED);
            return;
        }
    }
}

VariantContext::Type VariantContext::typeOfBiallelicVariant(const Allele & ref, const Allele & allele) {
    if ( ref.isSymbolic() )
        cerr << "Unexpected error: encountered a record with a symbolic reference allele" << endl;
    
    if ( allele.isSymbolic() )
        return SYMBOLIC;
    
    if ( ref.length() == allele.length() ) {
        if ( allele.length() == 1 )
            return SNP;
        else
            return MNP;
    }
    
    // Important note: previously we were checking that one allele is the prefix of the other.  However, that's not an
    // appropriate check as can be seen from the following example:
    // REF = CTTA and ALT = C,CT,CA
    // This should be assigned the INDEL type but was being marked as a MIXED type because of the prefix check.
    // In truth, it should be absolutely impossible to return a MIXED type from this method because it simply
    // performs a pairwise comparison of a single alternate allele against the reference allele (whereas the MIXED type
    // is reserved for cases of multiple alternate alleles of different types).  Therefore, if we've reached this point
    // in the code (so we're not a SNP, MNP, or symbolic allele), we absolutely must be an INDEL.
    return INDEL;
    
    // old incorrect logic:
    // if (oneIsPrefixOfOther(ref, allele))
    //     return Type.INDEL;
    // else
    //     return Type.MIXED;
}
#if 0

public String toString() {
    return String.format("[VC %s @ %s of type=%s alleles=%s attr=%s GT=%s",
                         getSource(), contig + ":" + (start - stop == 0 ? start : start + "-" + stop), this.getType(),
                         ParsingUtils.sortList(this.getAlleles()),
                         ParsingUtils.sortedString(this.getAttributes()),
                         this.getGenotypes());
}

// protected basic manipulation routines
private static List<Allele> makeAlleles(Collection<Allele> alleles) {
    final List<Allele> alleleList = new ArrayList<Allele>(alleles.size());
    
    boolean sawRef = false;
    for ( final Allele a : alleles ) {
        for ( final Allele b : alleleList ) {
            if ( a.equals(b, true) )
                throw new IllegalArgumentException("Duplicate allele added to VariantContext: " + a);
        }
        
        // deal with the case where the first allele isn't the reference
        if ( a.isReference() ) {
            if ( sawRef )
                throw new IllegalArgumentException("Alleles for a VariantContext must contain at most one reference allele: " + alleles);
            alleleList.add(0, a);
            sawRef = true;
        }
        else
            alleleList.add(a);
    }
    
    if ( alleleList.isEmpty() )
        throw new IllegalArgumentException("Cannot create a VariantContext with an empty allele list");
    
    if ( alleleList.get(0).isNonReference() )
        throw new IllegalArgumentException("Alleles for a VariantContext must contain at least one reference allele: " + alleles);
    
    return alleleList;
}

// ---------------------------------------------------------------------------------------------------------
//
// tribble integration routines -- not for public consumption
//
// ---------------------------------------------------------------------------------------------------------
#endif
    string VariantContext::getChr() {
    return contig;
}

int VariantContext::getStart() {
    return (int)start;
}

int VariantContext::getEnd() {
    return (int)stop;
}
#if 0

public boolean hasSymbolicAlleles() {
    for (final Allele a: getAlleles()) {
        if (a.isSymbolic()) {
            return true;
        }
    }
    return false;
}

public Allele getAltAlleleWithHighestAlleleCount() {
    // optimization: for bi-allelic sites, just return the 1only alt allele
    if ( isBiallelic() )
        return getAlternateAllele(0);
    
    Allele best = null;
    int maxAC1 = 0;
    for ( Allele a : getAlternateAlleles() ) {
        final int ac = getCalledChrCount(a);
        if ( ac >= maxAC1 ) {
            maxAC1 = ac;
            best = a;
        }
        
    }
    return best;
}

public int[] getGLIndecesOfAlternateAllele(Allele targetAllele) {
    
    int index = 1;
    for ( Allele allele : getAlternateAlleles() ) {
        if ( allele.equals(targetAllele) )
            break;
        index++;
    }
    
    return GenotypeLikelihoods.getPLIndecesOfAlleles(0, index);
}
#endif
