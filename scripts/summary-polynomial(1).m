if assigned batch then SetQuitOnError(true); else SetDebugOnError(true); end if;
AttachSpec("ASLoc.spec");
SetColumns(0);
SetAssertions(1);
SetDebugOnError(true);
SetVerbose("IHecke", 1);

procedure Usage(missingArg)
    printf "Error: %o\n", missingArg;
    print "";
    print "Usage example:";
    print "  magma -b type:=B5 prime:=2 saveDir:=saves summary-polynomial-KL.m > b5-2-polynomial-KL.txt";
    print "";
    print "Arguments:";
    print "  type             Required. Cartan type, e.g. A5, B5, C5, D5, E6.";
    print "  prime            Required. Prime characteristic p.";
    print "  targetLength     Optional. Only compute elements up to this length. Default: all.";
    print "  saveDir          Optional. Directory for saved .pcan files. Default: saves.";
    print "  maxWitnesses     Optional. Number of first changed elements to print. Default: 20.";
    print "  maxExamples      Optional. Number of largest correction examples to print in each table. Default: 20.";
    print "  maxPolyExamples  Optional. Number of polynomial frequency rows to print. Default: 20.";
    print "  printChanged     Optional. If true, print all changed elements. Default: false.";
    print "  printAllCorrectionStats Optional. If true, print one CSV row for every changed element. Default: false.";
    print "  printAllPolynomialStats Optional. If true, print one CSV row for every changed element with polynomial statistics. Default: false.";
    print "  chatty           Optional. ASLoc verbosity. Default: 1.";
    print "  profileName      Optional. Turn on profiling.";
    print "  stay             Optional. Stay around, do not quit.";
    quit;
end procedure;

// -----------------------------------------------------------------------------
// Arguments
// -----------------------------------------------------------------------------

if not assigned type then Usage("Missing argument: type"); end if;
cartanMat := CartanMatrix(type);
DynkinDiagram(cartanMat);
print "";

if not assigned prime then Usage("Missing argument: prime"); end if;
char := StringToInteger(prime);
if not IsPrime(char) then
    Usage("Invalid argument: prime should be a prime number.");
end if;

if not assigned saveDir then
    saveDir := "saves";
end if;

targetLength := assigned targetLength select StringToInteger(targetLength) else -1;
maxWitnesses := assigned maxWitnesses select StringToInteger(maxWitnesses) else 20;
maxExamples := assigned maxExamples select StringToInteger(maxExamples) else 20;
maxPolyExamples := assigned maxPolyExamples select StringToInteger(maxPolyExamples) else 20;
printChanged := assigned printChanged select (printChanged eq "true" or printChanged eq "True" or printChanged eq "1") else false;
printAllCorrectionStats := assigned printAllCorrectionStats select (printAllCorrectionStats eq "true" or printAllCorrectionStats eq "True" or printAllCorrectionStats eq "1") else false;
printAllPolynomialStats := assigned printAllPolynomialStats select (printAllPolynomialStats eq "true" or printAllPolynomialStats eq "True" or printAllPolynomialStats eq "1") else false;
chatty := assigned chatty select StringToInteger(chatty) else 1;
SetVerbose("ASLoc", chatty);

if assigned profileName then
    SetProfile(true);
    procedure WriteProfile()
        name := Sprintf("PCanPolynomialKLBasisSummary-Profile-%o-%o-%o", type, char, profileName);
        printf "Writing profile to %o.html\n", name;
        ProfileHTMLOutput(name);
    end procedure;
else
    procedure WriteProfile()
        printf "Profiling was not enabled.\n";
    end procedure;
end if;

// -----------------------------------------------------------------------------
// Formatting and small helpers
// -----------------------------------------------------------------------------

function FmtElt(w)
    return #w eq 0 select "id" else &cat[IntegerToString(i) : i in Eltseq(w)];
end function;

function FmtPoly(p)
    return Sprintf("%o", p);
end function;

function IsTruthyString(s)
    return s eq "true" or s eq "True" or s eq "TRUE" or s eq "1" or s eq "yes" or s eq "Yes";
end function;

function isPositiveSelfDual(g)
    LPoly<v> := BaseRing(Parent(g));
    supp, coeffs := Support(g);
    for p in coeffs do
        pcoeffs := Coefficients(p);
        if not forall{x : x in pcoeffs | x ge 0} then
            return false;
        end if;
        if p ne Evaluate(p, v^-1) then
            return false;
        end if;
    end for;
    return true;
end function;

function intersectLaurentPolys(f, g)
    LPoly := Parent(f);
    return &+[LPoly
        | Min(Coefficient(f, i), Coefficient(g, i))
        : i in [Min(Valuation(f), Valuation(g)) .. Max(Degree(f), Degree(g))]
    ];
end function;

function intersectBasisSummands(C, g, h)
    error if not isPositiveSelfDual(g), g, "is not positive self-dual";
    error if not isPositiveSelfDual(h), h, "is not positive self-dual";
    commonSupport := Support(g) meet Support(h);
    return &+[C | C.w * intersectLaurentPolys(Coefficient(g, w), Coefficient(h, w)) : w in commonSupport];
end function;

procedure InsertBasisElement(~pC, C, aut, w, pcan, eltsToCalculate)
    error if not isPositiveSelfDual(pcan), "not positive self dual";

    printf "Completed %o/%o: pC(%o) = %o\n", #pC, #eltsToCalculate, FmtElt(w), pcan;
    SetBasisElement(~pC, w, pcan);

    if w ne w^-1 and not IsDefined(pC, w^-1) then
        SetBasisElement(~pC, w^-1, &+[C | C.(y^-1) * Coefficient(pcan, y) : y in Support(pcan)]);
        printf "   Implied p-canonical basis element for inverse %o\n", FmtElt(w^-1);
    end if;

    if w ne aut(w) and not IsDefined(pC, aut(w)) then
        SetBasisElement(~pC, aut(w), &+[C | C.(aut(y)) * Coefficient(pcan, y) : y in Support(pcan)]);
        printf "   Implied p-canonical basis for diagram automorphism %o\n", FmtElt(aut(w));
    end if;

    if w ne aut(w^-1) and not IsDefined(pC, aut(w^-1)) then
        SetBasisElement(~pC, aut(w^-1), &+[C | C.(aut(y^-1)) * Coefficient(pcan, y) : y in Support(pcan)]);
        printf "   Implied p-canonical basis for diagram automorphism and inverse %o\n", FmtElt(aut(w^-1));
    end if;
end procedure;

procedure WriteBasis(pC, C, directory, type, p : complete:=false)
    LPoly<v> := BaseRing(pC);
    object := SerialiseBasis(pC);

    filename := Sprintf("%o/%o-%o.pcan", directory, type, p);
    if not complete then
        filename cat:= ".partial";
    end if;
    filenameTmp := filename cat ".tmp";

    fd := Open(filenameTmp, "w");
    fprintf fd, "%m", object;
    delete fd;
    System(Sprintf("mv \"%o\" \"%o\"", filenameTmp, filename));
    printf "%o p-canonical basis elements saved to %o\n", #pC, filename;
end procedure;

function ReadBasis(HAlg, directory, type, p)
    filename := Sprintf("%o/%o-%o.pcan", directory, type, p);
    ok, fd := OpenTest(filename, "rb");
    if not ok then
        filename := Sprintf("%o/%o-%o.pcan.partial", directory, type, p);
        ok, fd := OpenTest(filename, "rb");
        if not ok then
            return false, false;
        end if;
    end if;

    LPoly<v> := BaseRing(HAlg);
    pC := DeserialiseBasis(HAlg, eval Read(fd));
    delete fd;

    printf "%o p-canonical basis elements loaded from %o\n", #pC, filename;
    return pC, true;
end function;

function DiagramAut(W)
    id := hom<W -> W | [W.s : s in [1 .. Rank(W)]]>;

    if CartanName(W) eq "A~2" then
        return hom<W -> W | [W.2, W.1, W.3]>;
    end if;

    if not IsFinite(W) then
        return id;
    end if;

    if CartanName(W)[1] eq "A" then
        return hom<W -> W | [W.(Rank(W) - s + 1) : s in [1 .. Rank(W)]]>;
    elif CartanName(W)[1] eq "D" and Rank(W) ge 4 then
        return hom<W -> W | [W.s : s in [1 .. Rank(W) - 2] cat [Rank(W), Rank(W) - 1]]>;
    elif CartanName(W) eq "E6" then
        return hom<W -> W | [W.s : s in [6, 2, 5, 4, 3, 1]]>;
    end if;

    return id;
end function;


function CompareElts(x, y)
    if #x ne #y then
        return #x - #y;
    elif Eltseq(x) lt Eltseq(y) then
        return -1;
    elif Eltseq(x) eq Eltseq(y) then
        return 0;
    else
        return 1;
    end if;
end function;

function PolyMassAtOne(f)
    if f eq 0 then
        return 0;
    end if;
    coeffs := Coefficients(f);
    return #coeffs eq 0 select 0 else &+coeffs;
end function;

function PolyMaxCoefficient(f)
    if f eq 0 then
        return 0;
    end if;
    coeffs := Coefficients(f);
    return #coeffs eq 0 select 0 else Max(coeffs);
end function;

function PolyMaxAbsDegree(f)
    if f eq 0 then
        return 0;
    end if;
    degs := [d : d in [Valuation(f)..Degree(f)] | Coefficient(f, d) ne 0];
    return #degs eq 0 select 0 else Max([Abs(d) : d in degs]);
end function;


function PolyAbsMass(f)
    if f eq 0 then
        return 0;
    end if;
    coeffs := Coefficients(f);
    return #coeffs eq 0 select 0 else &+[Abs(c) : c in coeffs];
end function;

function PolyMonomialCount(f)
    if f eq 0 then
        return 0;
    end if;
    degs := [d : d in [Valuation(f)..Degree(f)] | Coefficient(f, d) ne 0];
    return #degs;
end function;

function PolyDegreeZeroCoeff(f)
    if f eq 0 then
        return 0;
    end if;
    return Coefficient(f, 0);
end function;

function PolyAbsDegreeZeroCoeff(f)
    return Abs(PolyDegreeZeroCoeff(f));
end function;

function PolyWeightedAbsDegreeMass(f)
    if f eq 0 then
        return 0;
    end if;
    return &+[Abs(d) * Abs(Coefficient(f, d)) : d in [Valuation(f)..Degree(f)] | Coefficient(f, d) ne 0];
end function;

function PolySecondDegreeMomentMass(f)
    if f eq 0 then
        return 0;
    end if;
    return &+[d*d * Abs(Coefficient(f, d)) : d in [Valuation(f)..Degree(f)] | Coefficient(f, d) ne 0];
end function;

function SafeRatio(num, den)
    RF := RealField(12);
    return den eq 0 select RF!0 else (RF!num)/(RF!den);
end function;

function CorrectionInCBasis(C, pC, w)
    // Measure the correction in the ordinary Kazhdan-Lusztig basis.
    // This is the positive convention:
    //     Delta_w = pC(w) - C(w) = sum_x f_{x,w}(v) C_x.
    return (C ! pC.w) - C.w;
end function;

function CorrectionStats(W, C, pC, w)
    // The correction is measured in the ordinary KL basis:
    //     Delta_w = pC(w) - C(w) = sum_x f_{x,w}(v) C(x).
    // This is the convention in which the correction coefficients should be
    // nonnegative Laurent polynomials.
    delta := CorrectionInCBasis(C, pC, w);
    supp := Sort(Setseq(Support(delta)));
    LPoly<v> := BaseRing(C);

    supportSize := #supp;
    if supportSize eq 0 then
        // Tuple format is documented below near the return statement.
        return <w, #w, 0, 0, 0, 0, 0, LPoly!0, 0, 0, 0, 0, 0, 0, 0, RealField(12)!0, RealField(12)!0, 0, 0>;
    end if;

    coeffs := [Coefficient(delta, x) : x in supp];

    masses := [PolyMassAtOne(f) : f in coeffs];
    absMasses := [PolyAbsMass(f) : f in coeffs];
    maxCoeffs := [PolyMaxCoefficient(f) : f in coeffs];
    maxDegrees := [PolyMaxAbsDegree(f) : f in coeffs];
    monomialCounts := [PolyMonomialCount(f) : f in coeffs];
    degreeZeroMasses := [PolyDegreeZeroCoeff(f) : f in coeffs];
    absDegreeZeroMasses := [PolyAbsDegreeZeroCoeff(f) : f in coeffs];
    weightedAbsDegreeMasses := [PolyWeightedAbsDegreeMass(f) : f in coeffs];
    secondDegreeMomentMasses := [PolySecondDegreeMomentMass(f) : f in coeffs];
    depths := [#w - #x : x in supp];

    gradedMassPolynomial := &+[LPoly | f : f in coeffs];
    signedMassAtOne := &+masses;
    absMassAtOne := &+absMasses;
    maxCoeff := Max(maxCoeffs);
    maxDegree := Max(maxDegrees);
    maxDepth := Max(depths);
    monomialComplexity := &+monomialCounts;
    degreeZeroMass := &+degreeZeroMasses;
    absDegreeZeroMass := &+absDegreeZeroMasses;
    nonzeroDegreeAbsMass := absMassAtOne - absDegreeZeroMass;
    weightedAbsDegreeMass := &+weightedAbsDegreeMasses;
    secondDegreeMomentMass := &+secondDegreeMomentMasses;
    degreeZeroRatio := SafeRatio(absDegreeZeroMass, absMassAtOne);
    avgAbsDegree := SafeRatio(weightedAbsDegreeMass, absMassAtOne);

    // Tuple format:
    //  1  w
    //  2  length(w)
    //  3  support_size in x
    //  4  signed_mass_at_one = sum_x f_x(1)        (old mass)
    //  5  max_coeff among coefficient polynomials  (old max coefficient)
    //  6  max_abs_degree                           (old degree spread)
    //  7  max_bruhat_depth
    //  8  graded_mass_polynomial = sum_x f_x(v)
    //  9  abs_mass_at_one = sum_{x,d} |[v^d]f_x|
    // 10  monomial_complexity = number of nonzero monomials across all f_x
    // 11  degree_zero_mass = sum_x [v^0]f_x
    // 12  abs_degree_zero_mass = sum_x |[v^0]f_x|
    // 13  nonzero_degree_abs_mass
    // 14  weighted_abs_degree_mass = sum_{x,d} |d| |[v^d]f_x|
    // 15  second_degree_moment_mass = sum_{x,d} d^2 |[v^d]f_x|
    // 16  degree_zero_ratio = abs_degree_zero_mass / abs_mass_at_one
    // 17  avg_abs_degree = weighted_abs_degree_mass / abs_mass_at_one
    // 18  max_abs_coeff = maximum absolute coefficient in all f_x
    // 19  number of distinct coefficient polynomials appearing in Delta_w
    maxAbsCoeff := Max([Max([Abs(c) : c in Coefficients(f)]) : f in coeffs]);
    distinctPolyStrings := {Sprintf("%o", f) : f in coeffs};

    return <w, #w, supportSize, signedMassAtOne, maxCoeff, maxDegree, maxDepth,
        gradedMassPolynomial, absMassAtOne, monomialComplexity, degreeZeroMass,
        absDegreeZeroMass, nonzeroDegreeAbsMass, weightedAbsDegreeMass,
        secondDegreeMomentMass, degreeZeroRatio, avgAbsDegree, maxAbsCoeff,
        #distinctPolyStrings>;
end function;

function StatAverage(stats, idx)
    RF := RealField(12);
    return #stats eq 0 select RF!0 else (RF!(&+[s[idx] : s in stats]))/(RF!#stats);
end function;

function StatMaximum(stats, idx)
    return #stats eq 0 select 0 else Max([s[idx] : s in stats]);
end function;

function StatCountMild(stats)
    // Old convention based on signed v=1 mass.
    return #[s : s in stats | s[4] eq 1];
end function;

function StatCountModerate(stats)
    return #[s : s in stats | s[4] ge 2 and s[4] le 5];
end function;

function StatCountWild(stats)
    return #[s : s in stats | s[4] gt 5];
end function;

function StatCountDegreeZeroOnly(stats)
    return #[s : s in stats | s[13] eq 0];
end function;

function StatCountGenuinelyGraded(stats)
    return #[s : s in stats | s[13] gt 0];
end function;

procedure IncrementStringCount(~A, key, inc)
    ok, val := IsDefined(A, key);
    A[key] := ok select val + inc else inc;
end procedure;

function StringCountRows(A)
    keys := Sort(Setseq(Keys(A)));
    rows := [<k, A[k]> : k in keys];
    Sort(~rows, func<a,b | (a[2] ne b[2]) select (b[2] - a[2]) else (a[1] lt b[1] select -1 else (a[1] eq b[1] select 0 else 1)) >);
    return rows;
end function;

procedure PrintTopCorrectionTable(title, stats, metricIndex, metricName, C, pC, maxExamples)
    print "";
    print title;
    printf "word,length,%o,support_size,signed_mass_at_one,abs_mass_at_one,monomial_complexity,max_coeff,max_abs_coeff,max_abs_degree,max_bruhat_depth,graded_mass_polynomial,pC(w)-C(w) in C-basis\n", metricName;

    Sort(~stats, func<a,b |
        (a[metricIndex] ne b[metricIndex]) select (b[metricIndex] - a[metricIndex])
        else CompareElts(a[1], b[1])
    >);

    nPrint := Min(maxExamples, #stats);
    for i in [1..nPrint] do
        s := stats[i];
        w := s[1];
        printf "%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o\n",
            FmtElt(w), s[2], s[metricIndex], s[3], s[4], s[9], s[10], s[5], s[18], s[6], s[7], s[8], CorrectionInCBasis(C, pC, w);
    end for;
end procedure;

procedure PrintPolynomialFrequencyTables(correctionStats, C, pC, maxPolyExamples)
    coeffPolyFreq := AssociativeArray();
    coeffPolyAbsMassFreq := AssociativeArray();
    gradedMassPolyFreq := AssociativeArray();

    for s in correctionStats do
        w := s[1];
        delta := CorrectionInCBasis(C, pC, w);
        supp := Sort(Setseq(Support(delta)));

        IncrementStringCount(~gradedMassPolyFreq, Sprintf("%o", s[8]), 1);

        for x in supp do
            f := Coefficient(delta, x);
            fstr := Sprintf("%o", f);
            IncrementStringCount(~coeffPolyFreq, fstr, 1);
            IncrementStringCount(~coeffPolyAbsMassFreq, fstr, PolyAbsMass(f));
        end for;
    end for;

    print "";
    print "COEFFICIENT POLYNOMIAL FREQUENCIES";
    print "polynomial,frequency,total_abs_mass";
    rows := StringCountRows(coeffPolyFreq);
    nPrint := Min(maxPolyExamples, #rows);
    for i in [1..nPrint] do
        r := rows[i];
        printf "%o,%o,%o\n", r[1], r[2], coeffPolyAbsMassFreq[r[1]];
    end for;

    print "";
    print "GRADED MASS POLYNOMIAL FREQUENCIES";
    print "graded_mass_polynomial,frequency";
    rows2 := StringCountRows(gradedMassPolyFreq);
    nPrint2 := Min(maxPolyExamples, #rows2);
    for i in [1..nPrint2] do
        r := rows2[i];
        printf "%o,%o\n", r[1], r[2];
    end for;
end procedure;

procedure PrintSummary(W, C, pC, type, char, targetLength, ifs, knownBySupports, torsionPrimes, formRecords, maxWitnesses, maxExamples, maxPolyExamples, printChanged, printAllCorrectionStats, printAllPolynomialStats)
    RF := RealField(12);
    allElts := Sort(Setseq(EnumerateCoxeterGroup(W : lengthBound := targetLength)));
    nTotal := #allElts;
    changed := [w : w in allElts | pC.w ne C.w];
    nChanged := #changed;
    pctChanged := (nTotal eq 0) select RF!0 else (RF!(100*nChanged))/(RF!nTotal);
    maxLen := Max([#w : w in allElts]);

    changedByLen := [#[w : w in changed | #w eq ell] : ell in [0..maxLen]];
    totalByLen := [#[w : w in allElts | #w eq ell] : ell in [0..maxLen]];

    firstLen := nChanged eq 0 select -1 else Min([#w : w in changed]);
    firstChanged := nChanged eq 0 select [] else [w : w in changed | #w eq firstLen];

    correctionStats := [CorrectionStats(W, C, pC, w) : w in changed];

    defectiveForms := [rec : rec in formRecords | rec[4] lt rec[5]];
    maxDefect := #defectiveForms eq 0 select 0 else Max([rec[5] - rec[4] : rec in defectiveForms]);

    avgSupport := StatAverage(correctionStats, 3);
    avgSignedMass := StatAverage(correctionStats, 4);
    avgAbsMass := StatAverage(correctionStats, 9);
    avgMonomialComplexity := StatAverage(correctionStats, 10);
    avgMaxCoeff := StatAverage(correctionStats, 5);
    avgMaxAbsCoeff := StatAverage(correctionStats, 18);
    avgMaxDegree := StatAverage(correctionStats, 6);
    avgMaxDepth := StatAverage(correctionStats, 7);
    avgDegreeZeroRatio := StatAverage(correctionStats, 16);
    avgAvgAbsDegree := StatAverage(correctionStats, 17);
    avgWeightedAbsDegreeMass := StatAverage(correctionStats, 14);

    maxSupport := StatMaximum(correctionStats, 3);
    maxSignedMass := StatMaximum(correctionStats, 4);
    maxAbsMass := StatMaximum(correctionStats, 9);
    maxMonomialComplexity := StatMaximum(correctionStats, 10);
    maxMaxCoeff := StatMaximum(correctionStats, 5);
    maxMaxAbsCoeff := StatMaximum(correctionStats, 18);
    maxMaxDegree := StatMaximum(correctionStats, 6);
    maxMaxDepth := StatMaximum(correctionStats, 7);
    maxWeightedAbsDegreeMass := StatMaximum(correctionStats, 14);
    maxAvgAbsDegree := StatMaximum(correctionStats, 17);

    mild := StatCountMild(correctionStats);
    moderate := StatCountModerate(correctionStats);
    wild := StatCountWild(correctionStats);
    degreeZeroOnly := StatCountDegreeZeroOnly(correctionStats);
    genuinelyGraded := StatCountGenuinelyGraded(correctionStats);

    print "";
    print "============================================================";
    print "PKL POLYNOMIAL SEVERITY SUMMARY";
    print "============================================================";
    printf "Type: %o\n", type;
    printf "Prime: %o\n", char;
    printf "Cartan name: %o\n", CartanName(W);
    printf "Rank: %o\n", Rank(W);
    printf "Target length: %o\n", targetLength eq -1 select "all" else IntegerToString(targetLength);
    printf "Elements considered: %o\n", nTotal;
    printf "Maximum length considered: %o\n", maxLen;
    printf "Different from KL: %o of %o\n", nChanged, nTotal;
    printf "Percentage different: %o%%\n", pctChanged;
    printf "Correction basis: ordinary KL basis, i.e. Delta_w = pC(w)-C(w) is expanded in the C-basis.\n";
    printf "Graded mass polynomial: M_w(v)=sum_x f_{x,w}(v), where Delta_w=sum_x f_{x,w}(v) C_x.\n";
    printf "Abs mass: sum of absolute values of all monomial coefficients in all f_{x,w}(v).\n";
    printf "Monomial complexity: total number of nonzero monomials across all correction polynomials f_{x,w}(v).\n";
    printf "Degree-zero ratio: abs degree-zero mass divided by abs mass.\n";
    printf "Mild/moderate/wild convention: signed_mass_at_one = 1 / 2..5 / >5.\n";
    printf "Elements certified by support-only shortcut: %o\n", knownBySupports;
    printf "Local intersection forms calculated during this run: %o\n", ifs;
    printf "Defective local intersection forms during this run: %o\n", #defectiveForms;
    printf "Maximum observed rank defect during this run: %o\n", maxDefect;
    printf "Torsion primes found during this run: %o\n", torsionPrimes;

    print "";
    print "OVERALL POLYNOMIAL SEVERITY";
    print "changed,avg_support_size,max_support_size,avg_signed_mass_at_one,max_signed_mass_at_one,avg_abs_mass_at_one,max_abs_mass_at_one,avg_monomial_complexity,max_monomial_complexity,avg_max_coeff,max_max_coeff,avg_max_abs_coeff,max_max_abs_coeff,avg_max_abs_degree,max_max_abs_degree,avg_max_bruhat_depth,max_max_bruhat_depth,avg_degree_zero_ratio,degree_zero_only,genuinely_graded,avg_avg_abs_degree,max_avg_abs_degree,avg_weighted_abs_degree_mass,max_weighted_abs_degree_mass,mild,moderate,wild";
    printf "%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o\n",
        nChanged, avgSupport, maxSupport, avgSignedMass, maxSignedMass, avgAbsMass, maxAbsMass,
        avgMonomialComplexity, maxMonomialComplexity, avgMaxCoeff, maxMaxCoeff,
        avgMaxAbsCoeff, maxMaxAbsCoeff, avgMaxDegree, maxMaxDegree, avgMaxDepth, maxMaxDepth,
        avgDegreeZeroRatio, degreeZeroOnly, genuinelyGraded, avgAvgAbsDegree, maxAvgAbsDegree,
        avgWeightedAbsDegreeMass, maxWeightedAbsDegreeMass, mild, moderate, wild;

    print "";
    print "LENGTH DISTRIBUTION";
    print "ell,total,changed,percentage_changed,avg_support_size,max_support_size,avg_signed_mass_at_one,max_signed_mass_at_one,avg_abs_mass_at_one,max_abs_mass_at_one,avg_monomial_complexity,max_monomial_complexity,avg_max_abs_degree,max_max_abs_degree,avg_max_bruhat_depth,max_max_bruhat_depth,avg_degree_zero_ratio,degree_zero_only,genuinely_graded,avg_avg_abs_degree,max_avg_abs_degree,avg_weighted_abs_degree_mass,max_weighted_abs_degree_mass,mild,moderate,wild";
    for ell in [0..maxLen] do
        statsAtLen := [s : s in correctionStats | s[2] eq ell];
        pct := totalByLen[ell+1] eq 0 select RF!0 else (RF!(100*changedByLen[ell+1]))/(RF!totalByLen[ell+1]);
        printf "%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o\n",
            ell,
            totalByLen[ell+1],
            changedByLen[ell+1],
            pct,
            StatAverage(statsAtLen, 3),
            StatMaximum(statsAtLen, 3),
            StatAverage(statsAtLen, 4),
            StatMaximum(statsAtLen, 4),
            StatAverage(statsAtLen, 9),
            StatMaximum(statsAtLen, 9),
            StatAverage(statsAtLen, 10),
            StatMaximum(statsAtLen, 10),
            StatAverage(statsAtLen, 6),
            StatMaximum(statsAtLen, 6),
            StatAverage(statsAtLen, 7),
            StatMaximum(statsAtLen, 7),
            StatAverage(statsAtLen, 16),
            StatCountDegreeZeroOnly(statsAtLen),
            StatCountGenuinelyGraded(statsAtLen),
            StatAverage(statsAtLen, 17),
            StatMaximum(statsAtLen, 17),
            StatAverage(statsAtLen, 14),
            StatMaximum(statsAtLen, 14),
            StatCountMild(statsAtLen),
            StatCountModerate(statsAtLen),
            StatCountWild(statsAtLen);
    end for;

    print "";
    print "FIRST WITNESSES";
    if nChanged eq 0 then
        print "No changed elements found.";
    else
        printf "First changed length: %o\n", firstLen;
        printf "Number changed at first changed length: %o\n", #firstChanged;
        Sort(~correctionStats, func<a,b | CompareElts(a[1], b[1]) >);
        nPrint := Min(maxWitnesses, #correctionStats);
        printf "Printing first %o changed elements by length/lex order.\n", nPrint;
        print "word,length,support_size,signed_mass_at_one,abs_mass_at_one,monomial_complexity,max_abs_degree,max_bruhat_depth,degree_zero_ratio,avg_abs_degree,graded_mass_polynomial,pC(w)-C(w) in C-basis";
        for i in [1..nPrint] do
            s := correctionStats[i];
            w := s[1];
            printf "%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o\n", FmtElt(w), s[2], s[3], s[4], s[9], s[10], s[6], s[7], s[16], s[17], s[8], CorrectionInCBasis(C, pC, w);
        end for;
    end if;

    if nChanged gt 0 then
        PrintTopCorrectionTable("LARGEST CORRECTIONS BY ABS MASS", correctionStats, 9, "abs_mass_at_one", C, pC, maxExamples);
        PrintTopCorrectionTable("LARGEST CORRECTIONS BY MONOMIAL COMPLEXITY", correctionStats, 10, "monomial_complexity", C, pC, maxExamples);
        PrintTopCorrectionTable("LARGEST CORRECTIONS BY SIGNED MASS", correctionStats, 4, "signed_mass_at_one", C, pC, maxExamples);
        PrintTopCorrectionTable("LARGEST CORRECTIONS BY DEGREE SPREAD", correctionStats, 6, "max_abs_degree", C, pC, maxExamples);
        PrintTopCorrectionTable("LARGEST CORRECTIONS BY WEIGHTED ABS DEGREE MASS", correctionStats, 14, "weighted_abs_degree_mass", C, pC, maxExamples);
        PrintTopCorrectionTable("DEEPEST CORRECTIONS IN BRUHAT LENGTH", correctionStats, 7, "max_bruhat_depth", C, pC, maxExamples);
        PrintPolynomialFrequencyTables(correctionStats, C, pC, maxPolyExamples);
    end if;

    print "";
    print "INTERSECTION FORM DEFECTS";
    print "word,x,degree,rank_mod_p,maxrank,defect,elementary_divisors";
    for rec in defectiveForms do
        w := rec[1];
        x := rec[2];
        deg := rec[3];
        rank := rec[4];
        maxrank := rec[5];
        elemDivs := rec[6];
        printf "%o,%o,%o,%o,%o,%o,%o\n", FmtElt(w), FmtElt(x), deg, rank, maxrank, maxrank-rank, elemDivs;
    end for;

    if printAllCorrectionStats or printAllPolynomialStats then
        print "";
        print "ALL POLYNOMIAL CORRECTION STATS";
        print "word,length,support_size,signed_mass_at_one,abs_mass_at_one,monomial_complexity,degree_zero_mass,abs_degree_zero_mass,nonzero_degree_abs_mass,weighted_abs_degree_mass,second_degree_moment_mass,degree_zero_ratio,avg_abs_degree,max_coeff,max_abs_coeff,max_abs_degree,max_bruhat_depth,distinct_coeff_polynomials,graded_mass_polynomial,pC(w)-C(w) in C-basis";
        Sort(~correctionStats, func<a,b | CompareElts(a[1], b[1]) >);
        for s in correctionStats do
            w := s[1];
            printf "%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o\n",
                FmtElt(w), s[2], s[3], s[4], s[9], s[10], s[11], s[12], s[13],
                s[14], s[15], s[16], s[17], s[5], s[18], s[6], s[7], s[19], s[8], CorrectionInCBasis(C, pC, w);
        end for;
    end if;

    if printChanged then
        print "";
        print "ALL CHANGED ELEMENTS";
        print "word,length,pC(w)-C(w) in C-basis";
        Sort(~changed, func<x,y | CompareElts(x, y) >);
        for w in changed do
            printf "%o,%o,%o\n", FmtElt(w), #w, CorrectionInCBasis(C, pC, w);
        end for;
    end if;

    print "============================================================";
    print "END PKL POLYNOMIAL SEVERITY SUMMARY";
    print "============================================================";
end procedure;

// -----------------------------------------------------------------------------
// Main calculation
// -----------------------------------------------------------------------------

W := CoxeterGroup(GrpFPCox, cartanMat);
B := BSParabolic(cartanMat, W, []);
HAlg, H, C := ShortcutIHeckeAlgebra(W);

eltsToCalculate := Sort(Setseq(EnumerateCoxeterGroup(W : lengthBound := targetLength)));

pC := CreateLiteralBasis(HAlg, "Canonical", "pC", Sprintf("%o-canonical basis of %o", char, type));
SetBasisElement(~pC, W.0, C.0);

if assigned saveDir then
    loaded, success := ReadBasis(HAlg, saveDir, type, char);
    if success then
        pC := loaded;
    end if;
end if;

aut := DiagramAut(W);

torsionPrimes := {Integers() |};
ifs := 0;
knownBySupports := 0;
formRecords := [];

lastSaveNum := #pC;
lastSaveTime := Realtime();

admissible_orders := {3} join (char ge 3 select {4} else {}) join (char ge 5 select {6} else {});
admissible_pairs := [{s, t} : s in [1 .. Rank(W)], t in [1 .. Rank(W)] | Order(W.s * W.t) in admissible_orders];

for i -> w in eltsToCalculate do
    if IsDefined(pC, w) then
        continue;
    end if;

    descElts :=
        [(C ! pC.(w * W.s)) * C.s : s in RightDescentSet(W, w)]
        cat
        [C.s * (C ! pC.(W.s * w)) : s in LeftDescentSet(W, w)];

    lrSupp := &meet[Support(h) : h in descElts];

    xdegsupp := &meet[
        {<x, deg> : deg in [0 .. Degree(Coefficient(h, x))], x in Support(h) | Coefficient(Coefficient(h, x), deg) ne 0}
        : h in descElts
    ];

    if #xdegsupp eq 1 then
        assert Rep(xdegsupp) eq <w, 0>;
        InsertBasisElement(~pC, C, aut, w, C.w, eltsToCalculate);
        knownBySupports +:= 1;
        continue;
    end if;

    degBoundTable := AssociativeArray();
    for h in descElts, x in Support(h), deg in [0 .. Degree(Coefficient(h, x))] do
        m_xwd := Coefficient(Coefficient(h, x), deg);
        ok, val := IsDefined(degBoundTable, <x, deg>);
        degBoundTable[<x, deg>] := ok select Min(val, m_xwd) else m_xwd;
    end for;

    knownCoeffs := AssociativeArray(W);
    for st in admissible_pairs do
        s, t := Explode(Setseq(st));
        wStar, ok := RightStar(W, w, s, t);
        if ok and #wStar lt #w then
            for y in Support(C ! pC.wStar) do
                yStar, ok := RightStar(W, y, s, t);
                if ok and yStar ne wStar and #yStar lt #w then
                    knownCoeffs[yStar] := Coefficient(C ! pC.wStar, y);
                end if;
            end for;
        end if;
    end for;
    for st in admissible_pairs do
        s, t := Explode(Setseq(st));
        wStar, ok := LeftStar(W, w, s, t);
        if ok and #wStar lt #w then
            for y in Support(C ! pC.wStar) do
                yStar, ok := LeftStar(W, y, s, t);
                if ok and yStar ne wStar and #yStar lt #w then
                    knownCoeffs[yStar] := Coefficient(C ! pC.wStar, y);
                end if;
            end for;
        end if;
    end for;

    if assigned saveDir and ((#pC gt lastSaveNum + 1000) or (Realtime(lastSaveTime) gt 60*5)) then
        WriteBasis(pC, C, saveDir, type, char);
        lastSaveNum := #pC;
        lastSaveTime := Realtime();
    end if;

    LPoly<v> := BaseRing(C);
    bigobj := &*[C.s : s in Eltseq(w)];

    printf "Calculating pcan for %o\n", C.w;

    xs := Reverse(Sort(Setseq(Support(bigobj))));
    for x in xs do
        if x notin lrSupp then
            bigobj -:= Coefficient(bigobj, x) * pC.x;
            continue;
        end if;

        if IsDefined(knownCoeffs, x) then
            bigobj := bigobj + (knownCoeffs[x] - Coefficient(bigobj, x)) * pC.x;
            continue;
        end if;

        error if LeftDescentSet(W, w) notsubset LeftDescentSet(W, x), "descent condition";
        error if RightDescentSet(W, w) notsubset RightDescentSet(W, x), "descent condition";

        lpoly := Coefficient(bigobj, x);
        for deg in [0 .. (lpoly eq 0 select -1 else Degree(lpoly))] do
            grade := deg eq 0 select 1 else v^deg + v^-deg;
            if Coefficient(lpoly, deg) eq 0 or x eq w then
                continue;
            end if;
            ok, bound := IsDefined(degBoundTable, <x, deg>);
            if not ok or bound eq 0 then
                bigobj := C ! bigobj - Coefficient(bigobj, x, deg) * grade * pC.x;
                continue;
            end if;

            maxrank := Coefficient(Coefficient(bigobj, x), deg);

            IF := LocalIntersectionForm(B, Eltseq(w), x, deg);
            ifs +:= 1;
            rank := Rank(ChangeRing(IF, GF(char)));

            elemDivs := ElementaryDivisors(IF);
            torsionPrimes join:= {Integers() | p : p in PrimeFactors(d), d in elemDivs};
            Append(~formRecords, <w, x, deg, rank, maxrank, elemDivs>);

            printf "   Degree %o intersection form of %o at %o has rank %o of maximum rank %o, elementary divisor set %o\n",
                deg, FmtElt(w), FmtElt(x), rank, maxrank, Seqset(elemDivs);

            bigobj := bigobj - rank * grade * pC.x;
        end for;
    end for;

    pcan := bigobj;
    InsertBasisElement(~pC, C, aut, w, pcan, eltsToCalculate);
end for;

if assigned saveDir and (#pC gt lastSaveNum) then
    WriteBasis(pC, C, saveDir, type, char : complete:=true);
    lastSaveNum := #pC;
    lastSaveTime := Realtime();
end if;

PrintSummary(W, C, pC, type, char, targetLength, ifs, knownBySupports, torsionPrimes, formRecords, maxWitnesses, maxExamples, maxPolyExamples, printChanged, printAllCorrectionStats, printAllPolynomialStats);

if assigned profileName then WriteProfile(); end if;
if not assigned stay then quit; end if;
