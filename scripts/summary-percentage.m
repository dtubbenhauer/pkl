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
    print "  magma -b type:=B5 prime:=2 saveDir:=saves percentage.m > b5-2.txt";
    print "";
    print "Arguments:";
    print "  type             Required. Cartan type, e.g. A5, B5, C5, D5, E6.";
    print "  prime            Required. Prime characteristic p.";
    print "  targetLength     Optional. Only compute elements up to this length. Default: all.";
    print "  saveDir          Optional. Directory for saved .pcan files. Default: saves.";
    print "  maxWitnesses     Optional. Number of first changed elements to print. Default: 20.";
    print "  printChanged     Optional. If true, print all changed elements. Default: false.";
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
printChanged := assigned printChanged select (printChanged eq "true" or printChanged eq "True" or printChanged eq "1") else false;
chatty := assigned chatty select StringToInteger(chatty) else 1;
SetVerbose("ASLoc", chatty);

if assigned profileName then
    SetProfile(true);
    procedure WriteProfile()
        name := Sprintf("PCanSummary-Profile-%o-%o-%o", type, char, profileName);
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

procedure PrintSummary(W, C, pC, type, char, targetLength, ifs, knownBySupports, torsionPrimes, formRecords, maxWitnesses, printChanged)
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

    // Defective local intersection forms, sorted by length of w, then degree.
    defectiveForms := [rec : rec in formRecords | rec[4] lt rec[5]];
    maxDefect := #defectiveForms eq 0 select 0 else Max([rec[5] - rec[4] : rec in defectiveForms]);

    print "";
    print "============================================================";
    print "PKL SUMMARY";
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
    printf "Elements certified by support-only shortcut: %o\n", knownBySupports;
    printf "Local intersection forms calculated: %o\n", ifs;
    printf "Defective local intersection forms: %o\n", #defectiveForms;
    printf "Maximum observed rank defect: %o\n", maxDefect;
    printf "Torsion primes found: %o\n", torsionPrimes;

    print "";
    print "LENGTH DISTRIBUTION";
    print "ell,total,changed,percentage_changed";
    for ell in [0..maxLen] do
        pct := totalByLen[ell+1] eq 0 select RF!0 else (RF!(100*changedByLen[ell+1]))/(RF!totalByLen[ell+1]);
        printf "%o,%o,%o,%o\n", ell, totalByLen[ell+1], changedByLen[ell+1], pct;
    end for;

    print "";
    print "FIRST WITNESSES";
    if nChanged eq 0 then
        print "No changed elements found.";
    else
        printf "First changed length: %o\n", firstLen;
        printf "Number changed at first changed length: %o\n", #firstChanged;
        Sort(~changed, func<x,y |
            (#x ne #y) select (#x - #y)
            else (Eltseq(x) lt Eltseq(y) select -1 else (Eltseq(x) eq Eltseq(y) select 0 else 1))
        >);
        nPrint := Min(maxWitnesses, #changed);
        printf "Printing first %o changed elements by length/lex order.\n", nPrint;
        print "word,length,pC(w)-C(w)";
        for i in [1..nPrint] do
            w := changed[i];
            printf "%o,%o,%o\n", FmtElt(w), #w, pC.w - C.w;
        end for;
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

    if printChanged then
        print "";
        print "ALL CHANGED ELEMENTS";
        print "word,length,pC(w)-C(w)";
        for w in changed do
            printf "%o,%o,%o\n", FmtElt(w), #w, pC.w - C.w;
        end for;
    end if;

    print "============================================================";
    print "END PKL SUMMARY";
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

PrintSummary(W, C, pC, type, char, targetLength, ifs, knownBySupports, torsionPrimes, formRecords, maxWitnesses, printChanged);

if assigned profileName then WriteProfile(); end if;
if not assigned stay then quit; end if;
