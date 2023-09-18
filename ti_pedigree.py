# A Python module for solving pedigree problems with a TI-Nspire CAS calculator


def pdsolve(pd, fmt="tbl"):
    """
    A random forest function for solving a pedigree given a nested list of its sub pedigrees (decision trees). The number of sub pedigrees must be equal to or larger than 2. The mode of inheritance is inferred entirely from the given/observed phenotypes (0s and/or 1s).
    The output of this function is either a list ("ls") or a table ("tbl"), and it should be specified using fmt.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    """

    pds = []
    if len(pd) < 2:
        return "Invalid input."
    for i in range(len(pd)):
        pds.append(pd[i].pdsolver(""))
    if not (fmt == "ls" or fmt == "tbl"):
        return "Invalid input."
    common = []
    autodom = 0
    autorec = 0
    sexdom = 0
    sexrec = 0
    mtdna = 0
    ychr = 0
    sumBoys = 0
    sumGirls = 0
    for i in range(len(pds[0])):
        sumBoys += sum(pd[0].boys)
        sumGirls += sum(pd[0].girls)
        if pds[0][i].find("autosomal dominant") != -1:
            common.append("autosomal dominant")
            if pds[0][i].find("+") != -1 and autodom == 0:
                autodom = "+"
            elif pds[0][i].find("+") != -1 and autodom == "-":
                autodom = "0"
            elif pds[0][i].find("-") != -1 and autodom == "+":
                autodom = "0"
        elif pds[0][i].find("autosomal recessive") != -1:
            common.append("autosomal recessive")
            if pds[0][i].find("+") != -1 and autorec == 0:
                autorec = "+"
            elif pds[0][i].find("+") != -1 and autorec == "-":
                autorec = 0
            elif pds[0][i].find("-") != -1 and autorec == "+":
                autorec = 0
        elif pds[0][i].find("sex-linked (X) dominant") != -1:
            common.append("sex-linked (X) dominant")
            if pds[0][i].find("++") != -1 and sexdom == 0:
                sexdom = "++"
            elif pds[0][i].find("+") != -1 and sexdom == 0:
                sexdom = "+"
            elif pds[0][i].find("+") != -1 and sexdom == "-":
                sexdom = 0
            elif pds[0][i].find("++") != -1 and sexdom == "-":
                sexdom = "+"
            elif pds[0][i][0] == "-" and sexdom == 0:
                sexdom = "-"
            elif pds[0][i][0] == "-" and sexdom == "+":
                sexdom = 0
            elif pds[0][i][0] == "-" and sexdom == "++":
                sexdom = "+"
        elif pds[0][i].find("sex-linked (X) recessive") != -1:
            common.append("sex-linked (X) recessive")
            if pds[0][i].find("+") != -1 and sexrec == 0:
                sexrec = "+"
            elif pds[0][i].find("+") != -1 and sexrec == "-":
                sexrec = 0
            elif pds[0][i][0] == "-" and sexrec == 0:
                sexrec = "-"
            elif pds[0][i][0] == "-" and sexrec == "+":
                sexrec = 0
        elif pds[0][i].find("mtDNA") != -1:
            common.append("mtDNA")
            if pds[0][i].find("+") != -1:
                if pds[0][i].find("++") != -1:
                    mtdna = "++"
                else:
                    mtdna = "+"
        elif pds[0][i].find("Y-chromosomal") != -1:
            common.append("Y-chromosomal")
            if pds[0][i].find("+") != -1:
                ychr = "+"
        # else:
        # 	return "Invalid input."

    for i in range(1, len(pds)):
        commonX = []
        sumBoys += sum(pd[i].boys)
        sumGirls += sum(pd[i].girls)
        for j in range(len(pds[i])):
            if pds[i][j].find("autosomal dominant") != -1:
                commonX.append("autosomal dominant")
                if pds[i][j].find("+") != -1 and autodom == 0:
                    autodom = "+"
                elif pds[i][j].find("+") != -1 and autodom == "-":
                    autodom = "0"
                elif pds[i][j].find("-") != -1 and autodom == "+":
                    autodom = "0"
            elif pds[i][j].find("autosomal recessive") != -1:
                commonX.append("autosomal recessive")
                if pds[i][j].find("+") != -1 and autorec == 0:
                    autorec = "+"
                elif pds[i][j].find("+") != -1 and autorec == "-":
                    autorec = 0
                elif pds[i][j].find("-") != -1 and autorec == "+":
                    autorec = 0
            elif pds[i][j].find("sex-linked (X) dominant") != -1:
                commonX.append("sex-linked (X) dominant")
                if pds[i][j].find("++") != -1 and sexdom == 0:
                    sexdom = "++"
                elif pds[i][j].find("+") != -1 and sexdom == 0:
                    sexdom = "+"
                elif pds[i][j].find("+") != -1 and sexdom == "-":
                    sexdom = 0
                elif pds[i][j].find("++") != -1 and sexdom == "-":
                    sexdom = "+"
                elif pds[i][j][0] == "-" and sexdom == 0:
                    sexdom = "-"
                elif pds[i][j][0] == "-" and sexdom == "+":
                    sexdom = 0
                elif pds[i][j][0] == "-" and sexdom == "++":
                    sexdom = "+"
            elif pds[i][j].find("sex-linked (X) recessive") != -1:
                commonX.append("sex-linked (X) recessive")
                if pds[i][j].find("+") != -1 and sexrec == 0:
                    sexrec = "+"
                elif pds[i][j].find("+") != -1 and sexrec == "-":
                    sexrec = 0
                elif pds[i][j][0] == "-" and sexrec == 0:
                    sexrec = "-"
                elif pds[i][j][0] == "-" and sexrec == "+":
                    sexrec = 0
            elif pds[i][j].find("mtDNA") != -1:
                commonX.append("mtDNA")
                if pds[i][j].find("+") != -1:
                    if pds[i][j].find("+") != -1:
                        if pds[i][j].find("++") != -1:
                            mtdna = "++"
                        else:
                            mtdna = "+"
            elif pds[i][j].find("Y-chromosomal") != -1:
                commonX.append("Y-chromosomal")
                if pds[i][j].find("+") != -1:
                    ychr = "+"
            # else:
            # 	return "Invalid input."

        common = list(set(common) & set(commonX))

    # Mode of inheritance (MofI).
    # According to the OMIM database autosomal recessive disorders are the most numerous, followed closely by autosomal dominant disorders. Sex-linked (X) recessive disorders are more numerous than sex-linked (X) dominant, which are rare. mtDNA and Y-chromosomal disorders are very rare. mtDNA disorders can follow any other inheritance pattern if the microdeletion site is in the nuclear mtDNA. The difference in observed numbers of autosomal and sex-linked disorders is mainly due to the different number of chromosomes that the chromosome types have (22 autosomal chromosomes and one sex chromosome in a haploid). The size of the chromosome does not seem to correlate with the number of disorders in that chromosome.
    # Prevalence of monogenic disorders in the general population (from Lynn B. Jorde et. al, Medical Genetics, 6th ed., Elsevier, 2020).
    # Autosomal dominant == 3/1000 to 9.5/1000
    # Autosomal recessive == 2/1000 to 2.5/1000
    # sex-linked (X) == 0.5/1000 to 2/1000
    # Y-chromosomal and mtDNA == rare
    MofI = []
    for i in range(len(common)):
        if common[i] == "sex-linked (X) dominant":
            if sexdom == "++" and sumGirls > sumBoys:
                MofI.append("++" + common[i])
            elif sexdom == "++":
                MofI.append("++" + common[i])
            # ?
            elif sexdom == "+" and sumGirls > sumBoys:
                MofI.append("++" + common[i])
            elif sexdom == "+":
                MofI.append("+" + common[i])
            elif sumGirls > sumBoys:
                MofI.append("+" + common[i])
            else:
                MofI.append(common[i])
        elif common[i] == "sex-linked (X) recessive":
            if sumBoys > sumGirls and sexrec == "+":
                MofI.append("++" + common[i])
            elif sexrec == "+":
                MofI.append("+" + common[i])
            elif sumBoys > sumGirls:
                MofI.append("+" + common[i])
            else:
                MofI.append(common[i])
        elif common[i] == "autosomal dominant":
            if autodom == "+":
                MofI.append("+" + common[i])
            else:
                MofI.append(common[i])
        elif common[i] == "autosomal recessive":
            if autorec == "+":
                MofI.append("+" + common[i])
            else:
                MofI.append(common[i])
        elif common[i] == "mtDNA":
            if mtdna == "+":
                MofI.append("+" + common[i])
            elif mtdna == "++":
                MofI.append("++" + common[i])
            else:
                MofI.append(common[i])
        elif common[i] == "Y-chromosomal":
            if ychr == "+" and sumBoys > 4 and sumGirls == 0:
                MofI.append("++" + common[i])
            elif ychr == "+":
                MofI.append("+" + common[i])
            elif sumBoys > 5 and sumGirls == 0:
                MofI.append("+" + common[i])
            elif sumBoys > 7 and sumGirls == 0:
                MofI.append("++" + common[i])
            else:
                MofI.append(common[i])
        else:
            MofI.append(common[i])

    for i in range(len(MofI)):
        if MofI[i][0] == "-":
            temp = MofI[i]
            MofI.pop(i)
            MofI.append(temp)

    for i in range(len(MofI)):
        if MofI[i][0] == "+":
            temp = MofI[i]
            MofI.pop(i)
            MofI.insert(0, temp)

    for i in range(len(MofI)):
        if MofI[i][0:2] == "++":
            temp = MofI[i]
            MofI.pop(i)
            MofI.insert(0, temp)

    if fmt == "ls":
        return MofI
    elif fmt == "tbl":
        print("Mode of inheritance:")
        for i in range(len(MofI)):
            print(MofI[i])


def prob(mofi, father=1, mother=1, hz="", sex="boy", status=1):
    """
    A function for estimating probabilities connected with a mode of inheritance of a pedigree
    -----------------------------------------------------------------------------------------------
    mofi: mode of inheritance:
                                                "ad": autosomal dominant
                                                "ar": autosomal recessive
                                                "sld": sex-linked (X) dominant
                                                "slr": sex-linked (X) recessive
                                                "mt": mitochondrial
                                                "y": Y-linked
    father: affected status of father (1 = affected, 0 = healthy)
    mother: affected status of mother (1 = affected, 0 = healthy)
    hz: heterozygosity of the parents:
                                                "": unknown
                                                "F hom": homozygous father
                                                "M hom": homozygous mother
                                                "F het": heterozygous father
                                                "M het": heterozygous mother
                                                "M&F hom": homozygous father and mother
                                                "F&M hom": homozygous father and mother
                                                "M&F het": heterozygous father and mother
                                                "F&M het": heterozygous father and mother
    sex: sex of the proband, boy or girl
    status: affected status of the proband (1 = affected, 0 = healthy)
    """

    # Check input.
    if not (
        mother == 1
        or mother == 0.5
        or mother == 0
        or father == 1
        or father == 0.5
        or father == 0
    ):
        return "Invalid input."
    if not (
        mofi == "ad"
        or mofi == "ar"
        or mofi == "sld"
        or mofi == "slr"
        or mofi == "mt"
        or mofi == "y"
    ):
        return "Invalid input."
    if not (sex == "boy" or sex == "girl"):
        return "Invalid input."
    if not (status == 1 or status == 0):
        return "Invalid input."
    if not (
        hz == ""
        or hz == "F hom"
        or hz == "M hom"
        or hz == "F het"
        or hz == "M het"
        or hz == "M&F hom"
        or hz == "F&M hom"
        or "M&F het"
        or "F&M het"
    ):
        return "Invalid input."
    prList = []
    if mofi == "ad":
        if father == 1 and mother == 1:
            # Homozygote parents, pr is 1.0.
            if hz == "F hom" or hz == "M hom":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            elif hz == "M het" or hz == "F het":
                if status == 1:
                    prList.append("autosomal dominant : 0.875")
                else:
                    prList.append("autosomal dominant : 0.125")
            # Het father and mother.
            elif hz == "M&F hom" or hz == "F&M hom":
                if status == 1:
                    prList.append("autosomal dominant : 0.75")
                else:
                    prList.append("autosomal dominant : 0.25")
            if hz == "M&F het" or hz == "F&M het":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            # Hom or het father/mother.
            else:
                if status == 1:
                    prList.append("autosomal dominant : 0.9375")
                else:
                    prList.append("autosomal dominant : 0.0625")
        elif father == 1 and mother == 0:
            # If father is hom then pr is 1.0.
            if hz == "F hom":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            # Het father.
            elif hz == "F het":
                if status == 1:
                    prList.append("autosomal dominant : 0.5")
                else:
                    prList.append("autosomal dominant : 0.5")
            # Hom or het father. Hom mother.
            elif hz == "" or hz == "M hom":
                if status == 1:
                    prList.append("autosomal dominant : 0.75")
                else:
                    prList.append("autosomal dominant : 0.25")
            elif hz == "M&F hom" or hz == "F&M hom":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            # hz == "M het" or hz == "M&F het" or hz == "F&M het"
            else:
                prList.append("autosomal dominant : N/A")
        elif father == 0 and mother == 1:
            # If mother is hom then pr is 1.0.
            if hz == "M hom":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            elif hz == "M het":
                if status == 1:
                    prList.append("autosomal dominant : 0.5")
                else:
                    prList.append("autosomal dominant : 0.5")
            # Hom or het mother. Hom father.
            elif hz == "" or hz == "F hom":
                if status == 1:
                    prList.append("autosomal dominant : 0.75")
                else:
                    prList.append("autosomal dominant : 0.25")
            elif hz == "M&F hom" or hz == "F&M hom":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            # hz == "M&F het" or hz == "F&M het"
            else:
                prList.append("autosomal dominant : N/A")
        elif father == 0.5 and mother == 1:
            # If mother is hom then pr is 1.
            if hz == "M hom":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            # Het mother.
            elif hz == "M het":
                if status == 1:
                    prList.append("autosomal dominant : 0.75")
                else:
                    prList.append("autosomal dominant : 0.25")
            # Hom or het mother.
            elif hz == "":
                if status == 1:
                    prList.append("autosomal dominant : 0.875")
                else:
                    prList.append("autosomal dominant : 0.125")
            elif hz == "M&F het" or hz == "F&M het":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            # hz == "F hom" or "M&F hom" or "F&M hom"
            else:
                prList.append("autosomal dominant : N/A")
        elif father == 1 and mother == 0.5:
            # F hom.
            if hz == "F hom":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            # Het mother.
            elif hz == "M het":
                if status == 1:
                    prList.append("autosomal dominant : 0.75")
                else:
                    prList.append("autosomal dominant : 0.25")
            # Hom or het father.
            elif hz == "":
                if status == 1:
                    prList.append("autosomal dominant : 0.875")
                else:
                    prList.append("autosomal dominant : 0.125")
            elif hz == "M&F het" or hz == "F&M het":
                if status == 1:
                    prList.append("autosomal dominant : 1.0")
                else:
                    prList.append("autosomal dominant : 0")
            # hz == "F hom" or "M&F hom" or "F&M hom"
            else:
                prList.append("autosomal dominant : N/A")
        elif father == 0.5 and mother == 0.5:
            if (
                hz == ""
                or hz == "M&F het"
                or hz == "F&M het"
                or hz == "M het"
                or hz == "F het"
            ):
                if status == 1:
                    prList.append("autosomal dominant : 0.75")
                else:
                    prList.append("autosomal dominant : 0.25")
            # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom"
            else:
                prList.append("autosomal dominant : N/A")
        elif father == 0 and mother == 0.5:
            if hz == "" or hz == "M het" or hz == "F hom":
                prList.append("autosomal dominant : 0.5")
            # hz == "M hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or hz == "M&F het" or hz == "F&M het"
            else:
                prList.append("autosomal dominant : N/A")
        elif father == 0.5 and mother == 0:
            if hz == "" or hz == "F het" or hz == "M hom":
                prList.append("autosomal dominant : 0.5")
            # hz == "M het" or hz == "F hom" or hz == "M&F hom" or hz == "F&M hom" or hz == "M&F het" or hz == "F&M het"
            else:
                prList.append("autosomal dominant : N/A")
        elif father == 0 and mother == 0:
            if (
                hz == ""
                or hz == "M hom"
                or hz == "F hom"
                or hz == "M&F hom"
                or hz == "F&M hom"
            ):
                if status == 1:
                    prList.append("autosomal dominant : 0")
                else:
                    prList.append("autosomal dominant : 1")
            # hz == "M het" or hz == "F het" or "M&F het" or "F&M het"
            else:
                prList.append("autosomal dominant : N/A")
    elif mofi == "ar":
        if father == 1 and mother == 1:
            if (
                hz == ""
                or hz == "M hom"
                or hz == "F hom"
                or hz == "M&F hom"
                or hz == "F&M hom"
            ):
                if status == 1:
                    prList.append("autosomal recessive : 1.0")
                else:
                    prList.append("autosomal recessive : 0")
            # hz == "M het" or hz == "F het" or hz == "M&F het" or hz == "F&M het"
            else:
                prList.append("autosomal recessive : N/A")
        elif father == 0.5 and mother == 1:
            if hz == "" or hz == "M hom" or hz == "F het":
                prList.append("autosomal recessive : 0.5")
            # hz == "M het" or hz == "F hom" or hz == "M&F het" or hz == "F&M het" or hz == "M&F hom" or hz == "F&M hom"
            else:
                prList.append("autosomal recessive : N/A")
        elif father == 1 and mother == 0.5:
            if hz == "" or hz == "M het" or hz == "F hom":
                prList.append("autosomal recessive : 0.5")
            # hz == "M hom" or hz == "F het" or hz == "M&F het" or hz == "F&M het" or hz == "M&F hom" or hz == "F&M hom"
            else:
                prList.append("autosomal recessive : N/A")
        elif father == 0.5 and mother == 0.5:
            if (
                hz == ""
                or hz == "M het"
                or hz == "F het"
                or hz == "M&F het"
                or hz == "F&M het"
            ):
                if status == 1:
                    prList.append("autosomal recessive : 0.25")
                else:
                    prList.append("autosomal recessive : 0.75")
            # hz == "M hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom"
            else:
                prList.append("autosomal recessive : N/A")
        elif father == 0.5 and mother == 0:
            # If mother is hom pr is 0.
            if hz == "M hom":
                if status == 1:
                    prList.append("autosomal recessive : 0")
                else:
                    prList.append("autosomal recessive : 1.0")
            # Het mother.
            elif hz == "M het" or hz == "M&F het" or hz == "F&M het":
                if status == 1:
                    prList.append("autosomal recessive : 0.25")
                else:
                    prList.append("autosomal recessive : 0.75")
            # Hom or het mother.
            elif hz == "" or hz == "F hom":
                if status == 1:
                    prList.append("autosomal recessive : 0.125")
                else:
                    prList.append("autosomal recessive : 0.875")
            # hz == "F hom" or hz == "M&F hom" or hz == "F&M hom"
            else:
                prList.append("autosomal recessive : N/A")
        elif father == 0 and mother == 0.5:
            # If father is hom pr is 0.
            if hz == "F hom":
                if status == 1:
                    prList.append("autosomal recessive : 0")
                else:
                    prList.append("autosomal recessive : 1.0")
            # Het father.
            elif hz == "F het" or hz == "M&F het" or hz == "F&M het":
                if status == 1:
                    prList.append("autosomal recessive : 0.25")
                else:
                    prList.append("autosomal recessive : 0.75")
            # Hom or het father.
            elif hz == "" or hz == "M hom":
                if status == 1:
                    prList.append("autosomal recessive : 0.125")
                else:
                    prList.append("autosomal recessive : 0.875")
            # hz == "M hom" or hz == "M&F hom" or hz == "F&M hom"
            else:
                prList.append("autosomal recessive : N/A")
        elif father == 1 and mother == 0:
            # If mother is hom pr is 0.
            if hz == "M hom":
                if status == 1:
                    prList.append("autosomal recessive : 0")
                else:
                    prList.append("autosomal recessive : 1.0")
            # Het mother.
            elif hz == "M het":
                if status == 1:
                    prList.append("autosomal recessive : 0.5")
                else:
                    prList.append("autosomal recessive : 0.5")
            elif hz == "M&F hom" or hz == "F&M hom":
                if status == 1:
                    prList.append("autosomal recessive : 0")
                else:
                    prList.append("autosomal recessive : 1.0")
            # Hom or het mother.
            elif hz == "" or hz == "F hom":
                if status == 1:
                    prList.append("autosomal recessive : 0.25")
                else:
                    prList.append("autosomal recessive : 0.75")
            # hz == "F het" or hz == "M&F het" or hz == "F&M het"
            else:
                prList.append("autosomal recessive : N/A")
        elif father == 0 and mother == 1:
            # If father is hom pr is 0.
            if hz == "F hom":
                if status == 1:
                    prList.append("autosomal recessive : 0")
                else:
                    prList.append("autosomal recessive : 1.0")
            # Het father.
            elif hz == "F het":
                if status == 1:
                    prList.append("autosomal recessive : 0.5")
                else:
                    prList.append("autosomal recessive : 0.5")
            elif hz == "M&F hom" or hz == "F&M hom":
                if status == 1:
                    prList.append("autosomal recessive : 0")
                else:
                    prList.append("autosomal recessive : 1.0")
            # Hom or het father.
            elif hz == "" or hz == "M hom":
                if status == 1:
                    prList.append("autosomal recessive : 0.25")
                else:
                    prList.append("autosomal recessive : 0.75")
            # hz == "M&F het" or hz == "F&M het"
            else:
                prList.append("autosomal recessive : N/A")
        elif father == 0 and mother == 0:
            # If mother or father is hom pr is 0.
            if hz == "F hom" or hz == "M hom":
                if status == 1:
                    prList.append("autosomal recessive : 0")
                else:
                    prList.append("autosomal recessive : 1.0")
            elif hz == "F het" or hz == "M het":
                if status == 1:
                    prList.append("autosomal recessive : 0.125")
                else:
                    prList.append("autosomal recessive : 0.875")
            # Het mother and father.
            elif hz == "M&F het" or hz == "F&M het":
                if status == 1:
                    prList.append("autosomal recessive : 0.25")
                else:
                    prList.append("autosomal recessive : 0.75")
            elif hz == "M&F hom" or hz == "F&M hom":
                if status == 1:
                    prList.append("autosomal recessive : 0")
                else:
                    prList.append("autosomal recessive : 1.0")
            # Hom or het father/mother.
            elif hz == "":
                if status == 1:
                    prList.append("autosomal recessive : 0.0625")
                else:
                    prList.append("autosomal recessive : 0.9375")
    elif mofi == "sld":
        if father == 1 and mother == 1:
            if sex == "boy":
                # If mother is hom pr is 1.0.
                if hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 1.0")
                    else:
                        prList.append("sex-linked (X) dominant : 0")
                elif hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 0.5")
                    else:
                        prList.append("sex-linked (X) dominant : 0.5")
                # Hom or het mother.
                elif hz == "":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 0.75")
                    else:
                        prList.append("sex-linked (X) dominant : 0.25")
                # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) dominant : N/A")
            elif sex == "girl":
                if hz == "" or hz == "M hom" or hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 1.0")
                    elif status == 0:
                        prList.append("sex-linked (X) dominant : 0")
                # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het":
                else:
                    prList.append("sex-linked (X) dominant : N/A")
        elif father == 1 and mother == 0.5:
            if sex == "boy":
                if hz == "" or hz == "M het":
                    prList.append("sex-linked (X) dominant : 0.5")
                else:
                    prList.append("sex-linked (X) dominant : N/A")
            elif sex == "girl":
                if hz == "" or hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 1.0")
                    elif status == 0:
                        prList.append("sex-linked (X) dominant : 0")
                # hz == "M hom" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) dominant : N/A")
        elif father == 1 and mother == 0:
            if sex == "boy":
                if status == 1:
                    prList.append("sex-linked (X) dominant : 0")
                elif status == 0:
                    prList.append("sex-linked (X) dominant : 1.0")
                elif (
                    hz == "M het"
                    or hz == "F hom"
                    or hz == "F het"
                    or hz == "M&F hom"
                    or hz == "F&M hom"
                    or "M&F het"
                    or "F&M het"
                ):
                    prList.append("sex-linked (X) dominant : N/A")
            elif sex == "girl":
                if status == 1:
                    prList.append("sex-linked (X) dominant : 1.0")
                elif status == 0:
                    prList.append("sex-linked (X) dominant : 0")
                elif (
                    hz == "M het"
                    or hz == "F hom"
                    or hz == "F het"
                    or hz == "M&F hom"
                    or hz == "F&M hom"
                    or "M&F het"
                    or "F&M het"
                ):
                    prList.append("sex-linked (X) dominant : N/A")
        elif father == 0 and mother == 1:
            if sex == "boy":
                # If mother is hom pr is 1.0.
                if hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 1.0")
                    elif status == 0:
                        prList.append("sex-linked (X) dominant : 0")
                # Het mother.
                elif hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 0.5")
                    elif status == 0:
                        prList.append("sex-linked (X) dominant : 0.5")
                # Hom or het mother.
                elif hz == "":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 0.75")
                    elif status == 0:
                        prList.append("sex-linked (X) dominant : 0.25")
                # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) dominant : N/A")
            elif sex == "girl":
                # If mother is hom pr is 1.0.
                if hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 1.0")
                    elif status == 0:
                        prList.append("sex-linked (X) dominant : 0")
                # Het mother.
                elif hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 0.5")
                    elif status == 0:
                        prList.append("sex-linked (X) dominant : 0.5")
                # Hom or het mother.
                elif hz == "":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 0.75")
                    else:
                        prList.append("sex-linked (X) dominant : 0.25")
                # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) dominant : N/A")
        elif father == 0 and mother == 0.5:
            if hz == "" or hz == "M het":
                if sex == "boy":
                    prList.append("sex-linked (X) dominant : 0.5")
                elif sex == "girl":
                    prList.append("sex-linked (X) dominant : 0.5")
            else:
                prList.append("sex-linked (X) dominant : N/A")
        elif father == 0 and mother == 0:
            if hz == "" or hz == "M hom":
                if sex == "boy":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 0")
                    else:
                        prList.append("sex-linked (X) dominant : 1.0")
                elif sex == "girl":
                    if status == 1:
                        prList.append("sex-linked (X) dominant : 0")
                    else:
                        prList.append("sex-linked (X) dominant : 1.0")
            # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
            else:
                prList.append("sex-linked (X) dominant : N/A")
        elif father == 0.5:
            prList.append("sex-linked (X) dominant : N/A")
    elif mofi == "slr":
        if father == 1 and mother == 1:
            if sex == "boy":
                if hz == "" or hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 1.0")
                    else:
                        prList.append("sex-linked (X) recessive : 0")
                # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
            elif sex == "girl":
                if hz == "" or hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 1.0")
                    else:
                        prList.append("sex-linked (X) recessive : 0")
                # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
        elif father == 1 and mother == 0.5:
            if hz == "" or hz == "M het":
                if sex == "boy":
                    prList.append("sex-linked (X) recessive : 0.5")
                elif sex == "girl":
                    prList.append("sex-linked (X) recessive : 0.5")
            # hz == "M hom" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
            else:
                prList.append("sex-linked (X) recessive : N/A")
        elif father == 1 and mother == 0:
            if sex == "boy":
                # If mother is hom then pr is 0.
                if hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0")
                    else:
                        prList.append("sex-linked (X) recessive : 1.0")
                # Het mother.
                elif hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0.5")
                    else:
                        prList.append("sex-linked (X) recessive : 0.5")
                # Hom or het mother.
                elif hz == "":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0.25")
                    else:
                        prList.append("sex-linked (X) recessive : 0.75")
                # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
            elif sex == "girl":
                # If mother is hom then pr is 0.
                if hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0")
                    else:
                        prList.append("sex-linked (X) recessive : 1.0")
                # Het mother.
                elif hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0.5")
                    else:
                        prList.append("sex-linked (X) recessive : 0.5")
                # Hom or het mother.
                elif hz == "":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0.25")
                    else:
                        prList.append("sex-linked (X) recessive : 0.75")
                # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
        elif father == 0 and mother == 1:
            if sex == "boy":
                if hz == "" or hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 1.0")
                    else:
                        prList.append("sex-linked (X) recessive : 0")
                # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
            elif sex == "girl":
                if hz == "" or hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0")
                    else:
                        prList.append("sex-linked (X) recessive : 1.0")
                # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
        elif father == 0 and mother == 0.5:
            if sex == "boy":
                if hz == "" or hz == "M het":
                    prList.append("sex-linked (X) recessive : 0.5")
                # hz == "M hom" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
            elif sex == "girl":
                if hz == "" or hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0")
                    else:
                        prList.append("sex-linked (X) recessive : 1.0")
                # hz == "M hom" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
        elif father == 0 and mother == 0:
            if sex == "boy":
                # If mother is hom then pr is 0.
                if hz == "M hom":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0")
                    else:
                        prList.append("sex-linked (X) recessive : 1.0")
                # Het mother.
                elif hz == "M het":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0.5")
                    else:
                        prList.append("sex-linked (X) recessive : 0.5")
                # Hom or het mother.
                elif hz == "":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0.25")
                    else:
                        prList.append("sex-linked (X) recessive : 0.75")
                # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
            elif sex == "girl":
                if hz == "" or hz == "M hom" or "M het":
                    if status == 1:
                        prList.append("sex-linked (X) recessive : 0")
                    else:
                        prList.append("sex-linked (X) recessive : 1.0")
                # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                else:
                    prList.append("sex-linked (X) recessive : N/A")
        elif father == 0.5:
            prList.append("sex-linked (X) recessive : N/A")
    elif mofi == "mt":
        # Maternal inheritance pattern.
        if hz == "":
            if mother == 1:
                if status == 1:
                    prList.append("mtDNA : 1.0")
                else:
                    prList.append("mtDNA : 0")
            elif mother == 0:
                prList.append("mtDNA : 0")
            else:
                prList.append("mtDNA : N/A")
        # Other inheritance patterns (autosomal dom./rec.,sex-linked (X) dom./rec., sporadic) connected with nuclear mtDNA disorders are not considered here. See the relevant sections: autosomal dom./rec. and sex-linked (X) dom./rec.
        # hz == "F hom" or hz == "M hom" or hz == "F het" or hz == "M het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het":
        else:
            prList.append("mtDNA : N/A")
    elif mofi == "y":
        # Paternal inheritance patttern.
        if hz == "":
            if father == 1 and mother == 0:
                if sex == "boy":
                    if status == 1:
                        prList.append("Y-chromosomal : 1.0")
                    else:
                        prList.append("Y-chromosomal : 0")
                elif sex == "girl":
                    if status == 1:
                        prList.append("Y-chromosomal : 0")
                    else:
                        prList.append("Y-chromosomal : 1.0")
            elif father == 0 and mother == 0:
                if sex == "boy":
                    if status == 1:
                        prList.append("Y-chromosomal : 0")
                    else:
                        prList.append("Y-chromosomal : 1")
                elif sex == "girl":
                    if status == 1:
                        prList.append("Y-chromosomal : 0")
                    else:
                        prList.append("Y-chromosomal : 1.0")
            else:
                prList.append("Y-chromosomal : N/A")
        # hz == "F hom" or hz == "M hom" or hz == "F het" or hz == "M het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het":
        else:
            prList.append("Y-chromosomal : N/A")

    if status == 1:
        st = "affected"
    else:
        st = "healthy"
    print("Probability of a " + str(sex) + " being " + st + ":")
    print(prList[0])


class pedigree:
    """
    ----------------------------------------------------------
    A Class for solving pedigree problems in medical genetics.
    ----------------------------------------------------------
    """

    def __init__(self, father, mother, boys, girls):
        self.father = father
        self.mother = mother
        self.boys = boys
        self.girls = girls

    def asdict(self):
        return {
            "father": self.father,
            "mother": self.mother,
            "boys": self.boys,
            "girls": self.girls,
        }

    def __str__(self):
        return (
            "father: "
            + str(self.father)
            + "\nmother: "
            + str(self.mother)
            + "\nboys: "
            + str(self.boys)
            + "\ngirls: "
            + str(self.girls)
        )

    def domrec(self):
        """
        A method that can be used to determine if a pedigree is dominant or recessive
        -----------------------------------------------------------------
        Dominant or recessive?
        Use [] for missing variables (only with boys/girls, father/mother cannot be missing).
        1. The parent's phenotypes have to be the same.
        2. At least some of the children's phenotypes should differ from their parents.
        """

        if self.mother != self.father or self.father == None or self.mother == None:
            return "Invalid input."

        if not (
            self.father == 0 or self.father == 1 or self.mother == 0 or self.mother == 1
        ):
            return "Invalid input."

        if len(self.boys) > 0:
            ZeroIndices = [i for i, x in enumerate(self.boys) if x == 0]
            OneIndices = [i for i, x in enumerate(self.boys) if x == 1]
            if len(ZeroIndices) + len(OneIndices) != len(self.boys):
                return "Invalid input."

        if len(self.girls) > 0:
            ZeroIndices = [i for i, x in enumerate(self.girls) if x == 0]
            OneIndices = [i for i, x in enumerate(self.girls) if x == 1]
            if len(ZeroIndices) + len(OneIndices) != len(self.girls):
                return "Invalid input."

        # Recessive.
        # Healthy parents have an affected child.
        if (
            self.father == 0
            and self.mother == 0
            and sum(self.girls) > 0
            and len(self.girls) != 0
        ):
            dr = "recessive"
        elif (
            self.father == 0
            and self.mother == 0
            and sum(self.boys) > 0
            and len(self.boys) != 0
        ):
            dr = "recessive"
        elif (
            self.father == 0
            and self.mother == 0
            and sum(self.boys) > 0
            and sum(self.girls) > 0
            and len(self.boys) != 0
            and len(self.girls) != 0
        ):
            dr = "recessive"
        # Parents' phenotype should be different from children's.
        elif (
            self.father == 0
            and self.mother == 0
            and sum(self.boys) == 0
            and sum(self.girls) == 0
            and len(self.boys) != 0
            and len(self.girls) != 0
        ):
            dr = "Invalid input."
            # "Undecided. The parents' phenotype should be different from their children's."
        elif (
            self.father == 0
            and self.mother == 0
            and sum(self.boys) == 0
            and sum(self.girls) == 0
            and len(self.boys) == 0
            and len(self.girls) == 0
        ):
            dr = "Invalid input."
            # "?"
        elif (
            self.father == 0
            and self.mother == 0
            and sum(self.boys) == 0
            and len(self.girls) == 0
        ):
            dr = "Invalid input."
            # "Undecided. The parents' phenotype should be different from their children's."
        elif (
            self.father == 0
            and self.mother == 0
            and len(self.boys) == 0
            and sum(self.girls) == 0
        ):
            dr = "Invalid input."
            # "Undecided. The parents' phenotype should be different from their children's."
        # Dominant.
        # Affected parents have a healthy child.
        elif (
            self.father == 1
            and self.mother == 1
            and sum(self.girls) < len(self.girls)
            and len(self.girls) != 0
        ):
            dr = "dominant"
        elif (
            self.father == 1
            and self.mother == 1
            and sum(self.boys) < len(self.boys)
            and len(self.boys) != 0
        ):
            dr = "dominant"
        elif (
            self.father == 1
            and self.mother == 1
            and (sum(self.boys) < len(self.boys) or sum(self.girls) < len(self.boys))
            and len(self.boys) != 0
            and len(self.girls) != 0
        ):
            dr = "dominant"
        # Parents' phenotype should be different from children's.
        elif (
            self.father == 1
            and self.mother == 1
            and sum(self.boys) == len(self.boys)
            and sum(self.girls) == 1
            and len(self.boys) != 0
            and len(self.girls) != 0
        ):
            dr = "Invalid input."
            # "Undecided. The parents' phenotype should be different from their children's."
        elif (
            self.father == 1
            and self.mother == 1
            and sum(self.boys) == len(self.boys)
            and len(self.girls) == 0
        ):
            dr = "Invalid input."
            # "Undecided. The parents' phenotype should be different from their children's."
        elif (
            self.father == 1
            and self.mother == 1
            and len(self.boys) == 0
            and sum(self.girls) == len(self.girls)
        ):
            dr = "Invalid input."
            # "Undecided. The parents' phenotype should be different from their children's."
        elif (
            self.father == 1
            and self.mother == 1
            and len(self.boys) == 0
            and len(self.girls) == 0
        ):
            dr = "Invalid input."
            # "?"
        else:
            return "Invalid input."

        return dr

    def sexLinked(self):
        # Analysis - sex-linked.

        # self.mother == 1
        if (
            self.mother == 1
            and sum(self.boys) == len(self.boys)
            and len(self.boys) != 0
        ):
            if len(self.girls) == 0:
                return ["autosomal recessive", "sex-linked (X) recessive"]
                # return "The disorder is autosomal or sex-linked (X) rec."
            elif sum(self.girls) == 0 and len(self.girls) != 0:
                return ["autosomal recessive", "sex-linked (X) recessive"]
                # return "The disorder is autosomal or sex-linked (X) rec."
            elif sum(self.girls) != 0 and len(self.girls) != 0:
                return ["autosomal recessive"]
                # return "The disorder is autosomal rec."
            else:
                return "? (sl)"
        # self.father == 1
        elif (
            self.father == 1
            and len(self.boys) == 0
            and sum(self.girls) == len(self.girls)
            and len(self.girls) != 0
        ):
            if len(self.girls) == 1:
                return [
                    "autosomal dominant",
                    "autosomal recessive",
                    "sex-linked (X) dominant",
                ]
                # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom."
            elif len(self.girls) == 2:
                return [
                    "+autosomal dominant",
                    "+sex-linked (X) dominant",
                    "autosomal recessive",
                ]
                # return "The disorder is more likely autosomal/sex-linked (X) dom. than autosomal rec."
            elif len(self.girls) == 3:
                return [
                    "+sex-linked (X) dominant",
                    "autosomal dominant",
                    "-autosomal recessive",
                ]
                # return "The disorder is more likely sex-linked (X) dom. than autosomal dom./rec."
            elif len(self.girls) > 3:
                return [
                    "++sex-linked (X) dominant",
                    "autosomal dominant",
                    "--autosomal recessive",
                ]
                # return "The disorder is more likely sex-linked (X) dom. than autosomal dom."
            else:
                return "? (sl)"
        elif (
            self.father == 1
            and sum(self.boys) == 0
            and sum(self.girls) == len(self.girls)
            and len(self.boys) != 0
            and len(self.girls) != 0
        ):
            if len(self.girls) <= 2:
                return ["sex-linked (X) dominant", "autosomal dominant"]
                # return "The disorder is sex-linked (X) or autosomal dom."
            elif len(self.girls) > 2:
                return ["+sex-linked (X) dominant", "autosomal dominant"]
                # return "The disorder is more likely sex-linked (X) dom. than autosomal dom."
            else:
                return "? (sl)"
        elif (
            self.father == 1
            and self.mother == 0
            and sum(self.boys) == len(self.boys)
            and len(self.boys) != 0
        ):
            if len(self.girls) == 0:
                if len(self.boys) < 3:
                    return [
                        "autosomal dominant",
                        "autosomal recessive (M het)",
                        "sex-linked (X) recessive (M het)",
                        "-Y-chromosomal",
                    ]
                    # return "The disorder could be autosomal dom., autosomal rec., sex-linked (X) rec. (heterozygous mother) or Y-chromosomal."
                elif len(self.boys) >= 3:
                    return [
                        "+Y-chromosomal",
                        "autosomal dominant",
                        "sex-linked (X) recessive (M het)",
                        "-autosomal recessive (M het)",
                    ]
                    # return "The disorder is likely Y-chromosomal, but it could be autosomal dom., autosomal or sex-linked (X) rec. (heterozygous mother)."
                else:
                    return "? (sl)"
            elif sum(self.girls) == 0 and len(self.girls) != 0:
                if len(self.boys) <= 2 and len(self.girls) <= 1:
                    return [
                        "autosomal dominant",
                        "autosomal recessive (M het)",
                        "sex-linked (X) recessive (M het)",
                        "-Y-chromosomal",
                    ]
                    # return "The disorder could be autosomal dom., autosomal or sex-linked (X) rec. (heterozygous mother) or Y-chromosomal."
                else:
                    return [
                        "+Y-chromosomal",
                        "autosomal dominant",
                        "autosomal recessive (M het)",
                        "-sex-linked (X) recessive (M het)",
                    ]
                    # return "The disorder is likely Y-chromosomal, but it could be autosomal dom., autosomal or sex-linked (X) rec. (heterozygous mother)."
                # No debugging breakpoint.
            elif sum(self.girls) != 0:
                return [
                    "autosomal dominant",
                    "autosomal recessive (M het)",
                    "sex-linked (X) recessive (M het)",
                ]
                # return "The disorder is autosomal dom., autosomal or sex-linked (X) rec. (heterozygous mother)."
            else:
                return "? (sl)"

    def ad(self):
        # Analysis - autosomal dominant.

        # self.father == 1 and self.mother == 0
        if self.father == 1 and self.mother == 0:
            if sum(self.girls) < len(self.girls) and len(self.girls) != 0:
                return ["autosomal dominant"]
                # return "The disorder is autosomal dom."
            elif sum(self.boys) > 0 and len(self.boys) != 0:
                return ["autosomal dominant"]
                # return "The disorder is autosomal dom."
            elif (
                sum(self.boys) == 0
                and sum(self.girls) == len(self.girls)
                and len(self.boys) != 0
                and len(self.girls) != 0
            ):
                # Sex-linked.
                return self.sexLinked()
            elif (
                sum(self.boys) == len(self.boys)
                and sum(self.girls) == len(self.girls)
                and len(self.boys) != 0
                and len(self.girls) != 0
            ):
                return ["autosomal dominant"]
                # return "The disorder is autosomal dom."
            elif (
                sum(self.boys) == 0
                and sum(self.girls) == len(self.girls)
                and len(self.boys) != 0
                and len(self.girls) != 0
            ):
                # Sex-linked.
                return self.sexLinked()
            else:
                return "? (ad)"
        # self.father == 0 and self.mother == 1
        elif self.father == 0 and self.mother == 1:
            if (
                sum(self.boys) <= len(self.boys)
                and sum(self.girls) <= len(self.girls)
                and len(self.boys) != 0
                and len(self.girls) != 0
            ):
                return ["autosomal dominant (M het)", "sex-linked (X) dominant (M het)"]
                # return "The disorder is autosomal or sex-linked (X) dom."
            elif (
                sum(self.boys) <= len(self.boys)
                and len(self.boys) != 0
                and len(self.girls) == 0
            ):
                return ["autosomal dominant (M het)", "sex-linked (X) dominant (M het)"]
                # return "The disorder is autosomal or sex-linked (X) dom."
            elif (
                len(self.boys) == 0
                and sum(self.girls) <= len(self.girls)
                and len(self.girls) != 0
            ):
                return ["autosomal dominant (M het)", "sex-linked (X) dominant (M het)"]
                # return "The disorder is autosomal or sex-linked (X) dom."
            else:
                return "? (ad)"
        # self.father == 1 and self.mother == 1
        elif self.father == 1 and self.mother == 1:
            if sum(self.girls) < len(self.girls) and len(self.girls) != 0:
                return ["autosomal dominant"]
                # return "The disorder is autosomal dom."
            elif (
                sum(self.boys) < len(self.boys)
                and len(self.boys) != 0
                and sum(self.girls) == len(self.girls)
                and len(self.girls) != 0
            ):
                return ["+sex-linked dominant", "autosomal dominant"]
                # return "The disorder is more likely sex-linked (X) than autosomal dom."
            elif (
                sum(self.boys) < len(self.boys)
                and len(self.boys) != 0
                and len(self.girls) == 0
            ):
                return ["+sex-linked dominant", "autosomal dominant"]
                # return "The disorder is more likely sex-linked (X) than autosomal dom."
            elif (
                sum(self.boys) == len(self.boys)
                and sum(self.girls) == len(self.girls)
                and len(self.boys) != 0
                and len(self.girls) != 0
            ):
                return ["+sex-linked (X) dominant", "autosomal dominant"]
                # return "The disorder is more likely sex-linked (X) than autosomal dom."
            elif sum(self.boys) == len(self.boys) and len(self.girls) == 0:
                return ["+autosomal dominant", "sex-linked (X) dominant"]
                # return "The disorder is more likely autosomal than sex-linked (X) dom."
            elif sum(self.girls) == len(self.girls) and len(self.boys) == 0:
                return ["sex-linked (X) dominant", "autosomal dominant"]
                # return "The disorder is sex-linked (X) or autosomal dom."
            else:
                return "? (ad)"

    def pdsolver(self, dr=""):
        """
        A decision tree for solving a pedigree. The mode of inheritance is inferred entirely from the observed/given phenotypes (0 or 1). The output of this function is a list of strings.
        ----------------------------------------------------------------------------------------------------------------------------------------------------
        dr: qualifier which determines if the mode of inheritance is dominant, recessive or unknown (dr = "dominant"/"recessive"/"")

        """

        # A decision tree for solving the mode of inheritance of pedigree problems.
        # Autosomal or sex-linked (X or Y)?
        # ...maybe mtDNA?

        # The parent's phenotypes should ideally be different from each other.

        # Use [] for missing variable(s) with boys and girls
        # Dominant or recessive?
        # dr = pedigree(father,mother,[boys],[girls]).domrec()
        # M = mother
        # F = father
        # het = heterozygous
        # hom = homozygous
        # +(+) = more likely
        # -(-) = less likely

        # Invalid input.
        if self.father == None or self.mother == None:
            return "Invalid input."

        if not (
            self.father == 0 or self.father == 1 or self.mother == 0 or self.mother == 1
        ):
            return "Invalid input."

        if (
            str(type(self.boys)) != "<class 'list'>"
            or str(type(self.girls)) != "<class 'list'>"
        ):
            return "Invalid input."

        if len(self.boys) > 0:
            ZeroIndices = [i for i, x in enumerate(self.boys) if x == 0]
            OneIndices = [i for i, x in enumerate(self.boys) if x == 1]
            if len(ZeroIndices) + len(OneIndices) != len(self.boys):
                return "Invalid input."

        if len(self.girls) > 0:
            ZeroIndices = [i for i, x in enumerate(self.girls) if x == 0]
            OneIndices = [i for i, x in enumerate(self.girls) if x == 1]
            if len(ZeroIndices) + len(OneIndices) != len(self.girls):
                return "Invalid input."

        # When you don't know if it's dominant or recessive, you have to guess.
        # Was it mitochondrial....?

        # 1
        if dr == "":
            # 2
            # self.father == 0 and self.mother == 1
            if self.father == 0 and self.mother == 1:
                # 3
                if sum(self.boys) != len(self.boys) and sum(self.girls) != len(
                    self.girls
                ):
                    # 4
                    # if sum(self.boys) <= len(self.boys) and sum(self.girls) <= len(self.girls):
                    # 5
                    if sum(self.boys) == 0 and sum(self.girls) == 0:
                        # indent 6?
                        if len(self.boys) == 0 and len(self.girls) == 0:
                            return [
                                "autosomal dominant",
                                "autosomal recessive",
                                "sex-linked (X) dominant",
                                "sex-linked(X) recessive",
                                "-mtDNA",
                            ]
                            # return "The disorder is autosomal or sex-linked (X) dom./rec. or mitochondrial."
                        # Redundancy ?
                        elif (
                            len(self.boys) + len(self.girls) == 1
                            and len(self.boys) == 0
                            and len(self.girls) != 0
                        ):
                            return [
                                "+autosomal recessive",
                                "autosomal dominant (M het)",
                                "sex-linked (X) dominant (M het)",
                            ]
                            # return "The disorder is autosomal recessive or autosomal/sex-linked (X) dom. (heterozygous mother)."
                        elif (
                            len(self.boys) + len(self.girls) == 1
                            and len(self.boys) != 0
                            and len(self.girls) == 0
                        ):
                            return [
                                "+autosomal recessive",
                                "autosomal dominant (M het)",
                                "sex-linked (X) dominant (M het)",
                            ]
                            # return "The disorder is autosomal recessive or autosomal/sex-linked dom. (heterozygous mother)."
                        elif (
                            len(self.boys) + len(self.girls) >= 2
                            and len(self.boys) != 0
                            and len(self.girls) != 0
                        ):
                            return [
                                "+autosomal recessive",
                                "autosomal dominant (M het)",
                                "sex-linked (X) dominant (M het)",
                            ]
                            # return "The disorder is autosomal recessive or autosomal/sex-linked dom. (heterozygous mother)."
                        elif (
                            len(self.boys) + len(self.girls) >= 2
                            and len(self.boys) == 0
                            and len(self.girls) != 0
                        ):
                            return [
                                "+autosomal recessive",
                                "autosomal dominant (M het)",
                                "sex-linked (X) dominant (M het)",
                            ]
                            # return "The disorder is autosomal recessive or autosomal/sex-linked dom. (heterozygous mother)."
                        elif (
                            len(self.boys) + len(self.girls) >= 2
                            and len(self.boys) != 0
                            and len(self.girls) == 0
                        ):
                            return [
                                "+autosomal recessive",
                                "autosomal dominant (M het)",
                                "sex-linked (X) dominant (M het)",
                            ]
                            # return "The disorder is autosomal recessive or autosomal/sex-linked dom. (heterozygous mother)."
                        else:
                            return "? (1)"
                        # indent 6?
                    # 5
                    elif (
                        sum(self.boys) + sum(self.girls) >= 1
                        and len(self.boys) != 0
                        or len(self.girls) != 0
                    ):
                        return [
                            "autosomal dominant (M het)",
                            "sex-linked (X) dominant (M het)",
                            "autosomal recessive (F het)",
                        ]
                        # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom."
                    else:
                        # Dominant.
                        # dr == "dominant"
                        return self.ad()
                        # 5
                    # 4
                # 3
                elif (
                    sum(self.boys) == len(self.boys)
                    and sum(self.girls) < len(self.girls)
                    and sum(self.girls) == 0
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return [
                        "+sex-linked (X) recessive",
                        "autosomal dominant",
                        "-sex-linked (X) dominant (M het)",
                        "-autosomal recessive (F het)",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                elif (
                    sum(self.boys) == len(self.boys)
                    and sum(self.girls) < len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return [
                        "+sex-linked (X) dominant",
                        "autosomal dominant",
                        "autosomal recessive",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom."
                elif (
                    sum(self.boys) < len(self.boys)
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.boys) != 0
                ):
                    return [
                        "+autosomal dominant",
                        "sex-linked (X) dominant",
                        "autosomal recessive",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom."
                # 3
                elif (
                    sum(self.boys) == len(self.boys)
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    # 4
                    if sum(self.boys) + sum(self.girls) < 3:
                        return [
                            "autosomal dominant",
                            "sex-linked (X) dominant",
                            "autosomal recessive (F het)",
                            "-mtDNA",
                        ]
                        # return "The disorder could be autosomal or sex-linked (X) dom, autosomal rec. (rec. if father is heterozygous) or mitochondrial."
                    elif sum(self.boys) + sum(self.girls) < 4:
                        return [
                            "+autosomal dominant",
                            "+sex-linked (X) dominant",
                            "autosomal recessive (F het)",
                            "mtDNA",
                        ]
                        # return "The disorder is more likely to be sex-linked/autosomal dom. than autosomal rec. or mitochondrial."
                    elif sum(self.boys) + sum(self.girls) >= 4:
                        return [
                            "++mtDNA",
                            "+autosomal dominant",
                            "sex-linked (X) dominant",
                            "--autosomal recessive (F het)",
                        ]
                        # return "The disorder could be inherited mitochondrially or it is sex-linked (X)/autosomal dom."
                    else:
                        return "? (2)"
                # 3
                elif (
                    sum(self.boys) == len(self.boys)
                    and len(self.boys) != 0
                    and len(self.girls) == 0
                ):
                    # 4
                    if sum(self.boys) < 4:
                        return [
                            "+autosomal dominant",
                            "+sex-linked (X) recessive (M hom)",
                            "autosomal recessive (F het)",
                            "sex-linked (X) dominant (M hom)",
                            "-mtDNA",
                        ]
                        # return "The disorder could be autosomal dom./rec. (rec. if father is heterozygous), sex-linked (X) dom./rec. or mitochondrial."
                    elif sum(self.boys) >= 4:
                        return [
                            "++mtDNA",
                            "autosomal dominant",
                            "-sex-linked (X) recessive (M hom)",
                            "-autosomal recessive (F het)",
                            "sex-linked (X) dominant (M hom)",
                        ]
                        # return "The disorder could be inherited mitochondrially (or it is autosomal dom.)."
                    else:
                        return "? (3)"
                    # 4
                # 3
                elif (
                    len(self.boys) == 0
                    and sum(self.girls) == len(self.girls)
                    and len(self.girls) != 0
                ):
                    # 4
                    if sum(self.girls) < 4:
                        return [
                            "+autosomal dominant",
                            "+sex-linked (X) dominant",
                            "autosomal recessive (F het)",
                            "-mtDNA",
                        ]
                        # return "The disorder could be autosomal dom./rec. (rec. if father is heterozygous), sex-linked (X) dominant or mitochondrial"
                    elif sum(self.girls) >= 4:
                        return ["+mtDNA", "autosomal dominant"]
                        # return "The disorder could be inherited mitochondrially (or it is autosomal dom.)."
                    else:
                        return "? (4)"
                    # 4
                elif (
                    sum(self.boys) == 0 and len(self.boys) != 0 and len(self.girls) == 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) dominant (M het)",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) recessive (M het)."
                elif (
                    sum(self.girls) == 0
                    and len(self.girls) != 0
                    and len(self.boys) == 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) dominant (M het)",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) recessive (M het)."
                # 3
                elif (
                    sum(self.boys) == len(self.boys)
                    and len(self.boys) != 0
                    and sum(self.girls) == len(self.girls)
                    and sum(self.boys) + sum(self.girls) >= 4
                    and len(self.girls) != 0
                ):
                    return ["+mtDNA", "-autosomal dominant"]
                    # return "The disorder could be inherited mitochondrially."
                else:
                    return "? (5)"
                # 3
            # self.father == 1 and self.mother == 0
            # 2
            elif self.father == 1 and self.mother == 0:
                # 3
                if sum(self.boys) == len(self.boys) and len(self.boys) != 0:
                    # Sex-linked.
                    return self.sexLinked()
                elif (
                    sum(self.boys) == len(self.boys)
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    # The mother has to be heterozygous for autosomal or sex-linked (X) rec.
                    # 4
                    if sum(self.boys) + sum(self.girls) < 3:
                        return [
                            "+autosomal dominant",
                            "autosomal recessive (M het)",
                            "sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder is autosomal dom. or autosomal/sex-linked (X) rec. (mother is heterozygous)."
                    elif sum(self.boys) + sum(self.girls) >= 3:
                        return [
                            "+autosomal dominant",
                            "-autosomal recessive (M het)",
                            "-sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder is more likely autosomal dom. than autosomal/sex-linked (X) rec. (mother is heterozygous)."
                    else:
                        return "? (6)"
                    # 4
                # 3
                elif (
                    sum(self.boys) == 0
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    # 4
                    if len(self.girls) <= 2:
                        return [
                            "+sex-linked (X) dominant",
                            "autosomal dominant (F het)",
                            "-autosomal recessive (M het)",
                            "-sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder is sex-linked (X) or autosomal dom. rather than autosomal or sex-linked (X) rec."
                    elif len(self.girls) > 2:
                        return [
                            "++sex-linked (X) dominant",
                            "autosomal dominant (F het)",
                            "-autosomal recessive (M het)",
                            "-sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder is sex-linked (X) or autosomal dom. rather than autosomal or sex-linked (X) rec."
                    else:
                        return "? (7)"
                    # 4
                # 3
                elif (
                    sum(self.girls) == len(self.girls)
                    and len(self.boys) == 0
                    and len(self.girls) != 0
                ):
                    # 4
                    if len(self.girls) == 1:
                        return [
                            "autosomal dominant",
                            "autosomal recessive",
                            "sex-linked (X) dominant",
                        ]
                        # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom."
                    elif len(self.girls) == 2:
                        return [
                            "+autosomal dominant",
                            "+sex-linked (X) dominant",
                            "autosomal recessive",
                        ]
                        # return "The disorder is more likely autosomal/sex-linked (X) dom. than autosomal rec."
                    elif len(self.girls) == 3:
                        return [
                            "+sex-linked (X) dominant",
                            "autosomal dominant",
                            "autosomal recessive",
                        ]
                        # return "The disorder is more likely sex-linked (X) dom. than autosomal dom./rec."
                    elif len(self.girls) > 3:
                        return [
                            "+sex-linked (X) dominant",
                            "autosomal dominant",
                            "--autosomal recessive",
                        ]
                        # return "The disorder is more likely sex-linked (X) dom. than autosomal dom."
                    else:
                        return "? (8)"
                    # 4
                elif (
                    sum(self.boys) == 0
                    and sum(self.girls) == 0
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) recessive",
                        "(--mtDNA)",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) rec."
                elif (
                    sum(self.boys) == 0
                    and sum(self.girls) == 0
                    and len(self.boys) == 0
                    and len(self.girls) == 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) dominant",
                        "sex-linked (X) recessive",
                        "(--Y-chromosomal)",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                elif (
                    sum(self.boys) == 0
                    and len(self.boys) == 0
                    and sum(self.girls) == 0
                    and len(self.girls) != 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) recessive",
                        "sex-linked (X) dominant",
                        "(--Y-chromosomal)",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                elif (
                    sum(self.girls) == 0
                    and len(self.girls) == 0
                    and sum(self.boys) == 0
                    and len(self.boys) != 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) recessive",
                        "(--mtDNA)",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) rec."
                # 3
                elif sum(self.boys) < len(self.boys) and sum(self.girls) < len(
                    self.girls
                ):
                    # 4
                    if len(self.boys) < 2 and len(self.girls) == 0:
                        return [
                            "autosomal dominant",
                            "autosomal recessive (M het)",
                            "sex-linked (X) recessive(M het)",
                        ]
                        # return "The disorder could be autosomal dom. or autosomal/sex-linked (X) rec. (heterozygous mother)."
                    elif len(self.girls) < 2 and len(self.boys) == 0:
                        return [
                            "autosomal dominant",
                            "autosomal recessive (M het)",
                            "sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder could be autosomal dom. or autosomal/sex-linked (X) rec. (heterozygous mother)."
                    elif len(self.boys) < 2 or len(self.girls) < 2:
                        return [
                            "autosomal dominant",
                            "autosomal recessive (M het)",
                            "sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder could be autosomal dom. or autosomal/sex-linked (X) rec. (heterozygous mother)."
                    # 4
                    elif len(self.boys) >= 2 or len(self.girls) >= 2:
                        # 5
                        if sum(self.boys) >= 2 or sum(self.girls) >= 2:
                            return [
                                "+autosomal dominant",
                                "autosomal recessive (M het)",
                                "sex-linked (X) recessive (M het)",
                            ]
                            # return "The disorder is more likely autosomal dom. than autosomal/sex-linked (X) rec. (heterozygous mother)."
                        else:
                            return [
                                "autosomal dominant",
                                "autosomal recessive (M het)",
                                "sex-linked (X) recessive (M het)",
                            ]
                            # return "The disorder could be autosomal dom. or autosomal/sex-linked (X) rec. (heterozygous mother."
                        # No debugging breakpoint.
                        # 5
                    # 4
                # 3
                elif sum(self.girls) < len(self.girls) and sum(self.boys) == len(
                    self.boys
                ):
                    # 4
                    if sum(self.girls) == 0 and len(self.girls) != 0:
                        # 5
                        if len(self.boys) == 0:
                            return [
                                "autosomal dominant (F het)",
                                "autosomal recessive",
                                "sex-linked (X) recessive",
                            ]
                            # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                        elif len(self.boys) == 1:
                            return [
                                "autosomal dominant (F het)",
                                "autosomal recessive",
                                "sex-linked (X) recessive",
                            ]
                            # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                        elif len(self.boys) > 1:
                            return [
                                "+autosomal dominant (F het)",
                                "autosomal recessive",
                                "sex-linked (X) recessive",
                            ]
                            # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                        else:
                            return "? (9)"
                        # 5
                    # 4
                    elif sum(self.girls) > 0:
                        # 5
                        # It is assumed that it is rare to be homozygous recessive with the disease (disorder,trait) allele.
                        return [
                            "+autosomal dominant (F het)",
                            "autosomal recessive (M het)",
                            "sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder is more likely autosomal dom. than autosomal or sex-linked (X) rec."
                    else:
                        return "? (10)"

                # 3
                elif sum(self.boys) < len(self.boys) and sum(self.girls) == len(
                    self.girls
                ):
                    # 4
                    if sum(self.boys) == 0 and len(self.boys) != 0:
                        # 5
                        if len(self.girls) == 0:
                            return [
                                "autosomal dominant (F het)",
                                "sex-linked (X) dominant (M hom)",
                                "autosomal recessive",
                                "sex-linked (X) recessive (M het)",
                            ]
                            # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                        elif len(self.girls) == 1:
                            return [
                                "sex-linked (X) dominant",
                                "autosomal dominant",
                                "autosomal recessive",
                                "sex-linked (X) recessive",
                            ]
                            # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                        elif len(self.girls) > 1:
                            return [
                                "+sex-linked (X) dominant",
                                "autosomal dominant",
                                "autosomal recessive",
                                "sex-linked (X) recessive",
                            ]
                            # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                        else:
                            return "? (11)"
                        # 5
                    # 4
                    elif sum(self.boys) > 0:
                        # 5
                        # It is assumed that it is rare to be homozygous recessive with the disease (disorder,trait) allele.
                        return [
                            "+autosomal dominant (F het)",
                            "autosomal recessive (M het)",
                            "sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder is more likely autosomal dom. than autosomal or sex-linked (X) rec."
                    else:
                        return "? (12)"

                    # 4:
                else:
                    return "? (13)"
                # 3
            # 2
            # self.father == 0 and self.mother == 0
            # 2
            elif self.father == 0 and self.mother == 0:
                # 3
                if sum(self.boys) > 0 and len(self.boys) != 0:
                    # 4
                    if len(self.girls) == 0:
                        return [
                            "+sex-linked (X) recessive (M het)",
                            "autosomal recessive (M&F het)",
                        ]
                        # return "The disorder is sex-linked (X) or autosomal rec."
                    elif sum(self.girls) == 0 and len(self.girls) != 0:
                        return [
                            "+sex-linked (X) recessive (M het)",
                            "autosomal recessive (M&F het)",
                        ]
                        # return "The disorder is more likely sex-linked (X) than autosomal rec."
                    elif sum(self.girls) > 0 and len(self.girls) != 0:
                        return ["autosomal recessive"]
                        # return "The disorder is autosomal rec."
                    else:
                        return "? (14)"
                # 3
                elif sum(self.girls) > 0:
                    return ["autosomal recessive"]
                    # return "The disorder is autosomal rec."
                elif sum(self.girls) == 0 and sum(self.boys) == 0:
                    # 4
                    # (--mtDNA)?
                    if len(self.boys) == 0:
                        return [
                            "autosomal dominant",
                            "autosomal recessive",
                            "sex-linked (X) dominant",
                            "sex-linked (X) recessive",
                            "(-Y-chromosomal)",
                        ]
                        # return "The disorder could be autosomal or sex-linked (X) dom./rec."
                    elif len(self.girls) == 0:
                        return [
                            "autosomal dominant",
                            "autosomal recessive",
                            "sex-linked (X) dominant",
                            "sex-linked (X) recessive",
                        ]
                        # return "The disorder could be autosomal or sex-linked (X) dom./rec."
                    else:
                        return [
                            "autosomal dominant",
                            "autosomal recessive",
                            "sex-linked (X) dominant",
                            "sex-linked (X) recessive",
                            "(-Y-chromosomal)",
                            "(--mtDNA)",
                        ]
                        # return "The disorder could be autosomal or sex-linked (X) dom./rec."
                    # No debugging breakpoint.
                    # 4
                # 3
                elif (
                    sum(self.boys) == 0
                    and sum(self.girls) == 0
                    and len(self.boys) == 0
                    and len(self.girls) == 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) dominant",
                        "sex-linked (X) recessive",
                        "(-Y-chromosomal)",
                        "(--mtDNA)",
                    ]
                    # return "The disorder could be autosomal dom./rec., sex-linked (X) dom./rec., mitochondrial or Y-chromosomal."
                else:
                    return "? (15)"
                # 3
            # 2
            # self.father == 1 and self.mother == 1
            # 2
            elif self.father == 1 and self.mother == 1:
                # 3
                if sum(self.girls) < len(self.girls) and len(self.girls) != 0:
                    return ["autosomal dominant"]
                    # return "The disorder is autosomal dom."
                elif (
                    sum(self.boys) < len(self.boys)
                    and len(self.boys) != 0
                    and sum(self.girls) == len(self.girls)
                    and len(self.girls) != 0
                ):
                    return ["+sex-linked (X) dominant", "autosomal dominant"]
                    # return "The disorder is sex-linked (X) or autosomal dom."
                elif (
                    sum(self.boys) < len(self.boys)
                    and len(self.boys) != 0
                    and len(self.girls) == 0
                ):
                    return ["+sex-linked (X) dominant", "autosomal dominant"]
                    # return "The disorder is sex-linked (X) or autosomal dom."
                elif (
                    sum(self.boys) == 0
                    and sum(self.girls) == 0
                    and len(self.boys) == 0
                    and len(self.girls) == 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) dominant",
                        "sex-linked (X) recessive",
                    ]
                    # return "The disorder could be autosomal dom./rec. or sex-linked (X) dom./rec."
                elif (
                    sum(self.boys) == len(self.boys)
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return [
                        "autosomal dominant",
                        "autosomal recessive",
                        "sex-linked (X) dominant",
                        "-sex-linked (X) recessive",
                    ]
                    # return "The disorder is autosomal or sex-linked (X) dom."
                elif sum(self.boys) == len(self.boys) and len(self.girls) == 0:
                    return ["+autosomal dominant", "sex-linked (X) dominant"]
                    # return "The disorder is autosomal or sex-linked (X) dom."
                elif sum(self.girls) == len(self.girls) and len(self.boys) == 0:
                    return ["sex-linked (X) dominant", "autosomal dominant"]
                    # return "The disorder is sex-linked (X) or autosomal dom."
                else:
                    return "? (16)"
                # 3
            # 2
        # 1
        # Dominant - autosomal or sex-linked
        # Autosomal

        if dr == "dominant":
            # self.father == 1 and self.mother == 0
            if self.father == 1 and self.mother == 0:
                if sum(self.girls) < len(self.girls) and len(self.girls) != 0:
                    return ["autosomal dominant"]
                    # return "The disorder is autosomal dom."
                elif sum(self.boys) > 0 and len(self.boys) != 0:
                    return ["autosomal dominant"]
                    # return "The disorder is autosomal dom."
                elif (
                    sum(self.boys) == 0
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    # Sex-linked.
                    return self.sexLinked()
                elif (
                    sum(self.boys) == len(self.boys)
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return ["autosomal dominant"]
                    # return "The disorder is autosomal dom."
                elif (
                    sum(self.boys) == 0
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    # Sex-linked.
                    return self.sexLinked()
                else:
                    return "? (15)"
            # self.father == 0 and self.mother == 1
            elif self.father == 0 and self.mother == 1:
                if (
                    sum(self.boys) <= len(self.boys)
                    and sum(self.girls) <= len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return [
                        "autosomal dominant (M het)",
                        "sex-linked (X) dominant (M het)",
                    ]
                    return "The disorder is autosomal or sex-linked (X) dom."
                elif (
                    sum(self.boys) <= len(self.boys)
                    and len(self.boys) != 0
                    and len(self.girls) == 0
                ):
                    return [
                        "autosomal dominant (M het)",
                        "sex-linked (X) dominant (M het)",
                    ]
                    # return "The disorder is autosomal or sex-linked (X) dom."
                elif (
                    len(self.boys) == 0
                    and sum(self.girls) <= len(self.girls)
                    and len(self.girls) != 0
                ):
                    return [
                        "autosomal dominant (M het)",
                        "sex-linked (X) dominant (M het)",
                    ]
                    # return "The disorder is autosomal or sex-linked (X) dom."
                else:
                    return "? (17)"
            # self.father == 1 and self.mother == 1
            elif self.father == 1 and self.mother == 1:
                if sum(self.girls) < len(self.girls) and len(self.girls) != 0:
                    return ["autosomal dominant"]
                    # return "The disorder is autosomal dom."
                elif (
                    sum(self.boys) < len(self.boys)
                    and len(self.boys) != 0
                    and sum(self.girls) == len(self.girls)
                    and len(self.girls) != 0
                ):
                    return ["+sex-linked dominant", "autosomal dominant"]
                    # return "The disorder is sex-linked (X) or autosomal dom."
                elif (
                    sum(self.boys) < len(self.boys)
                    and len(self.boys) != 0
                    and len(self.girls) == 0
                ):
                    return ["+sex-linked dominant", "autosomal dominant"]
                    # return "The disorder is likely sex-linked (X) autosomal dom."
                elif (
                    sum(self.boys) == len(self.boys)
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return ["+sex-linked (X) dominant", "autosomal dominant"]
                    return "The disorder is sex-linked (X) or autosomal dom."
                elif sum(self.boys) == len(self.boys) and len(self.girls) == 0:
                    return ["+autosomal dominant", "sex-linked (X) dominant"]
                    # return "The disorder is autosomal or sex-linked (X) dom."
                elif sum(self.girls) == len(self.girls) and len(self.boys) == 0:
                    return ["sex-linked (X) dominant", "autosomal dominant"]
                    # return "The disorder is sex-linked (X) or autosomal dom."
                else:
                    return "? (18)"

        # Recessive - autosomal or sex-linked
        # Autosomal

        # 1
        elif dr == "recessive":
            # self.father == 0 and self.mother == 1
            # 2
            if self.father == 0 and self.mother == 1:
                # 3
                if (
                    sum(self.boys) == 0
                    and sum(self.girls) == 0
                    and len(self.boys) == 0
                    and len(self.girls) == 0
                ):
                    return ["sex-linked (X) recessive", "autosomal recessive"]
                    # return "The disorder could be autosomal or sex-linked (X) rec."
                elif sum(self.boys) < len(self.boys) and len(self.boys) != 0:
                    return ["autosomal recessive"]
                    # return "The disorder is autosomal rec."
                elif sum(self.girls) != 0 and len(self.girls) != 0:
                    return ["autosomal recessive"]
                    # return "The disorder is autosomal rec."
                # 3
                elif len(self.girls) == 0:
                    # 4
                    if sum(self.boys) == len(self.boys):
                        # 5
                        if len(self.boys) < 2:
                            return ["autosomal recessive", "sex-linked (X) recessive"]
                            # return "The disorder could be autosomal or sex-linked (X) rec."
                        elif len(self.boys) >= 2:
                            return ["+sex-linked (X) recessive", "autosomal recessive"]
                            # return "The disorder is more likely sex-linked (X) than autosomal rec."
                        else:
                            return "? (19)"
                        # 5
                    # 4
                # 3
                elif len(self.boys) == 0:
                    # 4
                    # if sum(self.girls) != 0
                    if sum(self.girls) != 0 and len(self.girls) != 0:
                        return ["autosomal recessive"]
                        # return "The disorder is autosomal rec."
                    # 4
                    else:
                        # 5
                        # sum(self.girls) == 1
                        if len(self.girls) == 1 and len(self.girls) != 0:
                            return ["sex-linked (X) recessive", "autosomal recessive"]
                            # return "The disorder could be autosomal or sex-linked (X) rec."
                        # sum(self.girls) > 1
                        elif len(self.girls) > 1 and len(self.girls) != 0:
                            return ["+sex-linked (X) recessive", "autosomal recessive"]
                            # return "The disorder is more likely sex-linked (X) than autosomal rec."
                        else:
                            return "? (20)"
                        # 5
                    # 4
                # 3
                elif sum(self.boys) == len(self.boys):
                    # 4
                    # if sum(self.girls) == 1
                    if len(self.girls) == 1:
                        return ["sex-linked (X) recessive", "autosomal recessive"]
                        # return "The disorder could be autosomal or sex-linked (X) rec."
                    # if sum(self.girls) > 1
                    elif len(self.girls) > 1:
                        return ["+sex-linked (X) recessive", "autosomal recessive"]
                        # return "The disorder is more likely sex-linked (X) than autosomal rec."
                    # 4
                    else:
                        return "? (21)"
                else:
                    return "? (22)"
                # 3
                # ?
            # self.father == 1 and self.mother == 0
            # 2
            elif self.father == 1 and self.mother == 0:
                # 3
                if (
                    sum(self.girls) == 0
                    and sum(self.boys) == 0
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return ["autosomal recessive", "sex-linked (X) recessive"]
                    # return "The disorder is autosomal or sex-linked (X) rec."
                elif sum(self.boys) == 0 and sum(self.girls) == len(self.girls):
                    return [
                        "autosomal recessive (M het)",
                        "sex-linked (X) recessive (M het)",
                    ]
                    # return "The disorder is autosomal or sex-linked (X) rec. (heterozygote mother)."
                elif sum(self.boys) == len(self.boys) and sum(self.girls) == 0:
                    return [
                        "autosomal recessive (M het)",
                        "sex-linked (X) recessive (M het)",
                    ]
                    # return "The disorder is autosomal or sex-linked (X) rec. (heterozygote mother)."
                elif (
                    sum(self.boys) == len(self.boys)
                    and sum(self.girls) == len(self.girls)
                    and len(self.boys) != 0
                    and len(self.girls) != 0
                ):
                    return [
                        "autosomal recessive (M het)",
                        "sex-linked (X) recessive (M het)",
                    ]
                    # return "The disorder is autosomal or sex-linked (X) rec. (heterozygote mother)."
                elif sum(self.boys) < len(self.boys) and sum(self.girls) < len(
                    self.girls
                ):
                    # 4
                    if (
                        sum(self.boys) == 0
                        and sum(self.girls) == 0
                        and len(self.boys) == 0
                        and len(self.girls) == 0
                    ):
                        return ["autosomal recessive", "sex-linked (X) recessive"]
                        # return "The disorder could be autosomal or sex-linked (X) rec."
                    elif (
                        sum(self.boys) == 0
                        and len(self.boys) != 0
                        and sum(self.girls) == 0
                        and len(self.girls) != 0
                    ):
                        return ["autosomal recessive", "sex-linked (X) recessive"]
                        # return "The disorder could be autosomal or sex-linked (X) rec."
                    elif (
                        sum(self.boys) > 0
                        and len(self.boys) != 0
                        and sum(self.girls) > 0
                        and len(self.girls) != 0
                    ):
                        if sum(self.boys) > sum(self.girls):
                            return [
                                "+sex-linked (X) recessive (M het)",
                                "autosomal recessive (M het)",
                            ]
                            # return "The disorder is autosomal or sex-linked (X) rec. (heterozygote mother)."
                        else:
                            return [
                                "autosomal recessive (M het)",
                                "sex-linked (X) recessive (M het)",
                            ]
                            # return "The disorder is autosomal or sex-linked (X) rec. (heterozygote mother)."
                    elif sum(self.boys) > 0 and len(self.boys) != 0:
                        return [
                            "autosomal recessive (M het)",
                            "sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder is autosomal or sex-linked (X) rec. (heterozygote mother)."
                    elif sum(self.girls) > 0 and len(self.girls) != 0:
                        return [
                            "autosomal recessive (M het)",
                            "sex-linked (X) recessive (M het)",
                        ]
                        # return "The disorder is autosomal or sex-linked (X) rec. (heterozygote mother)."
                    else:
                        return "? (23)"
                    # 4
                # 3
                elif (
                    sum(self.boys) < len(self.boys)
                    and sum(self.boys) > 0
                    and sum(self.girls) == len(self.girls)
                ):
                    return [
                        "autosomal recessive (M het)",
                        "sex-linked (X) recessive (M het)",
                    ]
                    # return "The disorder is autosomal or sex-linked (X) rec. (heterozygote mother)."
                else:
                    return "? (24)"
                # 3
            # 2
            # self.father == 0 and self.mother == 0:
            elif self.father == 0 and self.mother == 0:
                # 3
                if sum(self.boys) > 0 and len(self.boys) != 0:
                    # 4
                    if len(self.girls) == 0:
                        return [
                            "+sex-linked (X) recessive (M het)",
                            "autosomal recessive (M&F het)",
                        ]
                        # return "The disorder is sex-linked (X) or autosomal rec."
                    elif sum(self.girls) == 0 and len(self.girls) != 0:
                        return [
                            "+sex-linked (X) recessive (M het)",
                            "autosomal recessive (M&F het)",
                        ]
                        # return "The disorder is more likely sex-linked (X) than autosomal rec."
                    elif sum(self.girls) > 0 and len(self.girls) != 0:
                        return ["autosomal recessive"]
                        # return "The disorder is autosomal rec."
                    else:
                        return "? (25)"
                    # 4
                # 3
                elif sum(self.girls) > 0:
                    return ["autosomal recessive"]
                    # return "The disorder is autosomal rec."
                elif sum(self.girls) == 0 and sum(self.boys) == 0:
                    # 4
                    if len(self.boys) == 0:
                        return ["+sex-linked (X) recessive", "autosomal recessive"]
                        # return "The disorder could be autosomal or sex-linked (X) rec."
                    elif len(self.girls) == 0:
                        return ["+autosomal recessive", "sex-linked (X) recessive"]
                        # return "The disorder could be autosomal or sex-linked (X) rec."
                    elif len(self.boys) != 0 and len(self.girls) != 0:
                        return ["autosomal recessive", "sex-linked (X) recessive"]
                        # return "The disorder could be autosomal or sex-linked (X) rec."
                    else:
                        return "? (26)"
                else:
                    return "? (27)"
            # 2
            else:
                # Sex-linked.
                return self.sexLinked()
            # 2
        # 1
        else:
            # Incorrect use of domrec() or invalid input.
            # If dr == "Invalid input.", dr == "Undecided. The parents' phenotype should be different from their children's.", dr == "?" or the input is invalid.
            return "Invalid input."

    def pr(self, dr="", hz="", sex="boy", status=1):
        """
        A method for tabulating average probabilities from a pedigree
        ------------------------------------------------------
        dr: qualifier which determines if the mode of inheritance is dominant, recessive or unknown (dr = "dominant"/"recessive"/"")
        hz: heterozygosity of the parents:
                        "": unknown
                        "F hom": homozygous father
                        "M hom": homozygous mother
                        "F het": heterozygous father
                        "M het": heterozygous mother
                        "M&F hom": homozygous father and mother
                        "F&M hom": homozygous father and mother
                        "M&F het": heterozygous father and mother
                        "F&M het": heterozygous father and mother
        sex: sex of the proband, boy or girl
        status: affected status of the proband (1 = affected, 0 = healthy/unaffected)
        """

        if not (dr == "" or dr == "dominant" or dr == "recessive"):
            return "Invalid input."

        pd = self.pdsolver(dr)

        if pd == "Invalid input.":
            return "Invalid input."
        if not (sex == "boy" or sex == "girl"):
            return "Invalid input."
        if not (status == 1 or status == 0):
            return "Invalid input."
        if not (
            hz == ""
            or hz == "F hom"
            or hz == "M hom"
            or hz == "F het"
            or hz == "M het"
            or hz == "M&F hom"
            or hz == "F&M hom"
            or "M&F het"
            or "F&M het"
        ):
            return "Invalid input."
        prList = []
        fatherHet = []
        motherHet = []
        for i in range(len(pd)):
            if (
                pd[i].find("M het") != -1
                or pd[i].find("M&F het") != -1
                or pd[i].find("F&M het") != -1
            ):
                motherHet.append(1)
            else:
                motherHet.append(0)
            if (
                pd[i].find("F het") != -1
                or pd[i].find("M&F het") != -1
                or pd[i].find("F&M het") != -1
            ):
                fatherHet.append(1)
            else:
                fatherHet.append(0)
        for i in range(len(pd)):
            if fatherHet[i] == 1:
                if hz != "":
                    father = 0.5
                else:
                    father = self.father
            else:
                father = self.father
            if motherHet[i] == 1:
                if hz != "":
                    mother = 0.5
                else:
                    mother = self.mother
            else:
                mother = self.mother
            if pd[i].find("autosomal dominant") != -1:
                # 1
                if father == 1 and mother == 1:
                    # Homozygote parents, pr is 1.0.
                    if hz == "F hom" or hz == "M hom":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    elif hz == "M het" or hz == "F het":
                        if status == 1:
                            prList.append("autosomal dominant : 0.875")
                        else:
                            prList.append("autosomal dominant : 0.125")
                    # Het father and mother.
                    elif hz == "M&F hom" or hz == "F&M hom":
                        if status == 1:
                            prList.append("autosomal dominant : 0.75")
                        else:
                            prList.append("autosomal dominant : 0.25")
                    if hz == "M&F het" or hz == "F&M het":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    # Hom or het father/mother.
                    else:
                        if status == 1:
                            prList.append("autosomal dominant : 0.9375")
                        else:
                            prList.append("autosomal dominant : 0.0625")
                # 2
                elif father == 1 and mother == 0:
                    # If father is hom then pr is 1.0.
                    if hz == "F hom":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    # Het father.
                    elif hz == "F het":
                        if status == 1:
                            prList.append("autosomal dominant : 0.5")
                        else:
                            prList.append("autosomal dominant : 0.5")
                    # Hom or het father. Hom mother.
                    elif hz == "" or hz == "M hom":
                        if status == 1:
                            prList.append("autosomal dominant : 0.75")
                        else:
                            prList.append("autosomal dominant : 0.25")
                    elif hz == "M&F hom" or hz == "F&M hom":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    # hz == "M het" or hz == "M&F het" or hz == "F&M het"
                    else:
                        prList.append("autosomal dominant : N/A")
                # 3
                elif father == 0 and mother == 1:
                    # If mother is hom then pr is 1.0.
                    if hz == "M hom":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    elif hz == "M het":
                        if status == 1:
                            prList.append("autosomal dominant : 0.5")
                        else:
                            prList.append("autosomal dominant : 0.5")
                    # Hom or het mother. Hom father.
                    elif hz == "" or hz == "F hom":
                        if status == 1:
                            prList.append("autosomal dominant : 0.75")
                        else:
                            prList.append("autosomal dominant : 0.25")
                    elif hz == "M&F hom" or hz == "F&M hom":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    # hz == "M&F het" or hz == "F&M het"
                    else:
                        prList.append("autosomal dominant : N/A")
                # 4
                elif father == 0.5 and mother == 1:
                    # If mother is hom then pr is 1.
                    if hz == "M hom":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    # Het mother.
                    elif hz == "M het":
                        if status == 1:
                            prList.append("autosomal dominant : 0.75")
                        else:
                            prList.append("autosomal dominant : 0.25")
                    # Hom or het mother.
                    elif hz == "":
                        if status == 1:
                            prList.append("autosomal dominant : 0.875")
                        else:
                            prList.append("autosomal dominant : 0.125")
                    elif hz == "M&F het" or hz == "F&M het":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    # hz == "F hom" or "M&F hom" or "F&M hom"
                    else:
                        prList.append("autosomal dominant : N/A")
                # 5
                elif father == 1 and mother == 0.5:
                    # If father is hom then pr is 1.
                    if hz == "F hom":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    # Het mother.
                    elif hz == "M het":
                        if status == 1:
                            prList.append("autosomal dominant : 0.75")
                        else:
                            prList.append("autosomal dominant : 0.25")
                    # Hom or het father.
                    elif hz == "":
                        if status == 1:
                            prList.append("autosomal dominant : 0.875")
                        else:
                            prList.append("autosomal dominant : 0.125")
                    elif hz == "M&F het" or hz == "F&M het":
                        if status == 1:
                            prList.append("autosomal dominant : 1.0")
                        else:
                            prList.append("autosomal dominant : 0")
                    # hz == "M hom" or "M&F hom" or "F&M hom"
                    else:
                        prList.append("autosomal dominant : N/A")
                # 6
                elif father == 0.5 and mother == 0.5:
                    if (
                        hz == ""
                        or hz == "M&F het"
                        or hz == "F&M het"
                        or hz == "M het"
                        or hz == "F het"
                    ):
                        if status == 1:
                            prList.append("autosomal dominant : 0.75")
                        else:
                            prList.append("autosomal dominant : 0.25")
                    # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom"
                    else:
                        prList.append("autosomal dominant : N/A")
                # 7
                elif father == 0 and mother == 0.5:
                    if hz == "" or hz == "M het" or hz == "F hom":
                        prList.append("autosomal dominant : 0.5")
                    # hz == "M hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or hz == "M&F het" or hz == "F&M het"
                    else:
                        prList.append("autosomal dominant : N/A")
                # 8
                elif father == 0.5 and mother == 0:
                    if hz == "" or hz == "F het" or hz == "M hom":
                        prList.append("autosomal dominant : 0.5")
                    # hz == "M het" or hz == "F hom" or hz == "M&F hom" or hz == "F&M hom" or hz == "M&F het" or hz == "F&M het"
                    else:
                        prList.append("autosomal dominant : N/A")
                # 9
                elif father == 0 and mother == 0:
                    if (
                        hz == ""
                        or hz == "M hom"
                        or hz == "F hom"
                        or hz == "M&F hom"
                        or hz == "F&M hom"
                    ):
                        if status == 1:
                            prList.append("autosomal dominant : 0")
                        else:
                            prList.append("autosomal dominant : 1")
                    # hz == "M het" or hz == "F het" or "M&F het" or "F&M het"
                    else:
                        prList.append("autosomal dominant : N/A")
            elif pd[i].find("autosomal recessive") != -1:
                # 1
                if father == 1 and mother == 1:
                    if (
                        hz == ""
                        or hz == "M hom"
                        or hz == "F hom"
                        or hz == "M&F hom"
                        or hz == "F&M hom"
                    ):
                        if status == 1:
                            prList.append("autosomal recessive : 1.0")
                        else:
                            prList.append("autosomal recessive : 0")
                    # hz == "M het" or hz == "F het" or hz == "M&F het" or hz == "F&M het"
                    else:
                        prList.append("autosomal recessive : N/A")
                # 2
                elif father == 0.5 and mother == 1:
                    if hz == "" or hz == "M hom" or hz == "F het":
                        prList.append("autosomal recessive : 0.5")
                    # hz == "M het" or hz == "F hom" or hz == "M&F het" or hz == "F&M het" or hz == "M&F hom" or hz == "F&M hom"
                    else:
                        prList.append("autosomal recessive : N/A")
                # 3
                elif father == 1 and mother == 0.5:
                    if hz == "" or hz == "M het" or hz == "F hom":
                        prList.append("autosomal recessive : 0.5")
                    # hz == "M hom" or hz == "F het" or hz == "M&F het" or hz == "F&M het" or hz == "M&F hom" or hz == "F&M hom"
                    else:
                        prList.append("autosomal recessive : N/A")
                # 4
                elif father == 0.5 and mother == 0.5:
                    if (
                        hz == ""
                        or hz == "M het"
                        or hz == "F het"
                        or hz == "M&F het"
                        or hz == "F&M het"
                    ):
                        if status == 1:
                            prList.append("autosomal recessive : 0.25")
                        else:
                            prList.append("autosomal recessive : 0.75")
                    # hz == "M hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom"
                    else:
                        prList.append("autosomal recessive : N/A")
                # 5
                elif father == 0.5 and mother == 0:
                    # If mother is hom pr is 0.
                    if hz == "M hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0")
                        else:
                            prList.append("autosomal recessive : 1.0")
                    # Het mother.
                    elif hz == "M het" or hz == "M&F het" or hz == "F&M het":
                        if status == 1:
                            prList.append("autosomal recessive : 0.25")
                        else:
                            prList.append("autosomal recessive : 0.75")
                    # Hom or het mother.
                    elif hz == "" or hz == "F hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0.125")
                        else:
                            prList.append("autosomal recessive : 0.875")
                    # hz == "F hom" or hz == "M&F hom" or hz == "F&M hom"
                    else:
                        prList.append("autosomal recessive : N/A")
                # 6
                elif father == 0 and mother == 0.5:
                    # If father is hom pr is 0.
                    if hz == "F hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0")
                        else:
                            prList.append("autosomal recessive : 1.0")
                    # Het father.
                    elif hz == "F het" or hz == "M&F het" or hz == "F&M het":
                        if status == 1:
                            prList.append("autosomal recessive : 0.25")
                        else:
                            prList.append("autosomal recessive : 0.75")
                    # Hom or het father.
                    elif hz == "" or hz == "M hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0.125")
                        else:
                            prList.append("autosomal recessive : 0.875")
                    # hz == "M hom" or hz == "M&F hom" or hz == "F&M hom"
                    else:
                        prList.append("autosomal recessive : N/A")
                # 7
                elif father == 1 and mother == 0:
                    # If mother is hom pr is 0.
                    if hz == "M hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0")
                        else:
                            prList.append("autosomal recessive : 1.0")
                    # Het mother.
                    elif hz == "M het":
                        if status == 1:
                            prList.append("autosomal recessive : 0.5")
                        else:
                            prList.append("autosomal recessive : 0.5")
                    elif hz == "M&F hom" or hz == "F&M hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0")
                        else:
                            prList.append("autosomal recessive : 1.0")
                    # Hom or het mother.
                    elif hz == "" or hz == "F hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0.25")
                        else:
                            prList.append("autosomal recessive : 0.75")
                    # hz == "F het" or hz == "M&F het" or hz == "F&M het"
                    else:
                        prList.append("autosomal recessive : N/A")
                # 8
                elif father == 0 and mother == 1:
                    # If father is hom pr is 0.
                    if hz == "F hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0")
                        else:
                            prList.append("autosomal recessive : 1.0")
                    # Het father.
                    elif hz == "F het":
                        if status == 1:
                            prList.append("autosomal recessive : 0.5")
                        else:
                            prList.append("autosomal recessive : 0.5")
                    elif hz == "M&F hom" or hz == "F&M hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0")
                        else:
                            prList.append("autosomal recessive : 1.0")
                    # Hom or het father.
                    elif hz == "" or hz == "M hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0.25")
                        else:
                            prList.append("autosomal recessive : 0.75")
                    # hz == "M&F het" or hz == "F&M het"
                    else:
                        prList.append("autosomal recessive : N/A")
                # 9
                elif father == 0 and mother == 0:
                    # If mother or father is hom pr is 0.
                    if hz == "F hom" or hz == "M hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0")
                        else:
                            prList.append("autosomal recessive : 1.0")
                    elif hz == "F het" or hz == "M het":
                        if status == 1:
                            prList.append("autosomal recessive : 0.125")
                        else:
                            prList.append("autosomal recessive : 0.875")
                    # Het mother and father.
                    elif hz == "M&F het" or hz == "F&M het":
                        if status == 1:
                            prList.append("autosomal recessive : 0.25")
                        else:
                            prList.append("autosomal recessive : 0.75")
                    elif hz == "M&F hom" or hz == "F&M hom":
                        if status == 1:
                            prList.append("autosomal recessive : 0")
                        else:
                            prList.append("autosomal recessive : 1.0")
                    # Hom or het father/mother.
                    elif hz == "":
                        if status == 1:
                            prList.append("autosomal recessive : 0.0625")
                        else:
                            prList.append("autosomal recessive : 0.9375")
            elif pd[i].find("sex-linked (X) dominant") != -1:
                # 1
                if father == 1 and mother == 1:
                    if sex == "boy":
                        # If mother is hom pr is 1.0.
                        if hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 1.0")
                            else:
                                prList.append("sex-linked (X) dominant : 0")
                        elif hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 0.5")
                            else:
                                prList.append("sex-linked (X) dominant : 0.5")
                        # Hom or het mother.
                        elif hz == "":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 0.75")
                            else:
                                prList.append("sex-linked (X) dominant : 0.25")
                        # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) dominant : N/A")
                    elif sex == "girl":
                        if hz == "" or hz == "M hom" or hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 1.0")
                            elif status == 0:
                                prList.append("sex-linked (X) dominant : 0")
                        # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het":
                        else:
                            prList.append("sex-linked (X) dominant : N/A")
                # 2
                elif father == 1 and mother == 0.5:
                    if sex == "boy":
                        if hz == "" or hz == "M het":
                            prList.append("sex-linked (X) dominant : 0.5")
                        else:
                            prList.append("sex-linked (X) dominant : N/A")
                    elif sex == "girl":
                        if hz == "" or hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 1.0")
                            elif status == 0:
                                prList.append("sex-linked (X) dominant : 0")
                        # hz == "M hom" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) dominant : N/A")
                # 3
                elif father == 1 and mother == 0:
                    if sex == "boy":
                        if status == 1:
                            prList.append("sex-linked (X) dominant : 0")
                        elif status == 0:
                            prList.append("sex-linked (X) dominant : 1.0")
                        elif (
                            hz == "M het"
                            or hz == "F hom"
                            or hz == "F het"
                            or hz == "M&F hom"
                            or hz == "F&M hom"
                            or "M&F het"
                            or "F&M het"
                        ):
                            prList.append("sex-linked (X) dominant : N/A")
                    elif sex == "girl":
                        if status == 1:
                            prList.append("sex-linked (X) dominant : 1.0")
                        elif status == 0:
                            prList.append("sex-linked (X) dominant : 0")
                        elif (
                            hz == "M het"
                            or hz == "F hom"
                            or hz == "F het"
                            or hz == "M&F hom"
                            or hz == "F&M hom"
                            or "M&F het"
                            or "F&M het"
                        ):
                            prList.append("sex-linked (X) dominant : N/A")
                # 4
                elif father == 0 and mother == 1:
                    if sex == "boy":
                        # If mother is hom pr is 1.0.
                        if hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 1.0")
                            elif status == 0:
                                prList.append("sex-linked (X) dominant : 0")
                        # Het mother.
                        elif hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 0.5")
                            elif status == 0:
                                prList.append("sex-linked (X) dominant : 0.5")
                        # Hom or het mother.
                        elif hz == "":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 0.75")
                            elif status == 0:
                                prList.append("sex-linked (X) dominant : 0.25")
                        # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) dominant : N/A")
                    elif sex == "girl":
                        # If mother is hom pr is 1.0.
                        if hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 1.0")
                            elif status == 0:
                                prList.append("sex-linked (X) dominant : 0")
                        # Het mother.
                        elif hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 0.5")
                            elif status == 0:
                                prList.append("sex-linked (X) dominant : 0.5")
                        # Hom or het mother.
                        elif hz == "":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 0.75")
                            else:
                                prList.append("sex-linked (X) dominant : 0.25")
                        # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) dominant : N/A")
                # 5
                elif father == 0 and mother == 0.5:
                    if hz == "" or hz == "M het":
                        if sex == "boy":
                            prList.append("sex-linked (X) dominant : 0.5")
                        elif sex == "girl":
                            prList.append("sex-linked (X) dominant : 0.5")
                    else:
                        prList.append("sex-linked (X) dominant : N/A")
                # 6
                elif father == 0 and mother == 0:
                    if hz == "" or hz == "M hom":
                        if sex == "boy":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 0")
                            else:
                                prList.append("sex-linked (X) dominant : 1.0")
                        elif sex == "girl":
                            if status == 1:
                                prList.append("sex-linked (X) dominant : 0")
                            else:
                                prList.append("sex-linked (X) dominant : 1.0")
                    # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                    else:
                        prList.append("sex-linked (X) dominant : N/A")
            elif pd[i].find("sex-linked (X) recessive") != -1:
                # 1
                if father == 1 and mother == 1:
                    if sex == "boy":
                        if hz == "" or hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 1.0")
                            else:
                                prList.append("sex-linked (X) recessive : 0")
                        # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                    elif sex == "girl":
                        if hz == "" or hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 1.0")
                            else:
                                prList.append("sex-linked (X) recessive : 0")
                        # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                # 2
                elif father == 1 and mother == 0.5:
                    if hz == "" or hz == "M het":
                        if sex == "boy":
                            prList.append("sex-linked (X) recessive : 0.5")
                        elif sex == "girl":
                            prList.append("sex-linked (X) recessive : 0.5")
                    # hz == "M hom" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                    else:
                        prList.append("sex-linked (X) recessive : N/A")
                # 3
                elif father == 1 and mother == 0:
                    if sex == "boy":
                        # If mother is hom then pr is 0.
                        if hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0")
                            else:
                                prList.append("sex-linked (X) recessive : 1.0")
                        # Het mother.
                        elif hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0.5")
                            else:
                                prList.append("sex-linked (X) recessive : 0.5")
                        # Hom or het mother.
                        elif hz == "":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0.25")
                            else:
                                prList.append("sex-linked (X) recessive : 0.75")
                        # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                    elif sex == "girl":
                        # If mother is hom then pr is 0.
                        if hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0")
                            else:
                                prList.append("sex-linked (X) recessive : 1.0")
                        # Het mother.
                        elif hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0.5")
                            else:
                                prList.append("sex-linked (X) recessive : 0.5")
                        # Hom or het mother.
                        elif hz == "":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0.25")
                            else:
                                prList.append("sex-linked (X) recessive : 0.75")
                        # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                # 4
                elif father == 0 and mother == 1:
                    if sex == "boy":
                        if hz == "" or hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 1.0")
                            else:
                                prList.append("sex-linked (X) recessive : 0")
                        # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                    elif sex == "girl":
                        if hz == "" or hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0")
                            else:
                                prList.append("sex-linked (X) recessive : 1.0")
                        # hz == "M het" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                # 5
                elif father == 0 and mother == 0.5:
                    if sex == "boy":
                        if hz == "" or hz == "M het":
                            prList.append("sex-linked (X) recessive : 0.5")
                        # hz == "M hom" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                    elif sex == "girl":
                        if hz == "" or hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0")
                            else:
                                prList.append("sex-linked (X) recessive : 1.0")
                        # hz == "M hom" or hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                # 6
                elif father == 0 and mother == 0:
                    if sex == "boy":
                        # If mother is hom then pr is 0.
                        if hz == "M hom":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0")
                            else:
                                prList.append("sex-linked (X) recessive : 1.0")
                        # Het mother.
                        elif hz == "M het":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0.5")
                            else:
                                prList.append("sex-linked (X) recessive : 0.5")
                        # Hom or het mother.
                        elif hz == "":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0.25")
                            else:
                                prList.append("sex-linked (X) recessive : 0.75")
                        # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
                    elif sex == "girl":
                        if hz == "" or hz == "M hom" or "M het":
                            if status == 1:
                                prList.append("sex-linked (X) recessive : 0")
                            else:
                                prList.append("sex-linked (X) recessive : 1.0")
                        # hz == "F hom" or hz == "F het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het"
                        else:
                            prList.append("sex-linked (X) recessive : N/A")
            elif pd[i].find("mtDNA") != -1:
                # Maternal inheritance pattern.
                if hz == "":
                    if mother == 1:
                        if status == 1:
                            prList.append("mtDNA : 1.0")
                        else:
                            prList.append("mtDNA : 0")
                    elif mother == 0:
                        if status == 1:
                            prList.append("mtDNA : 0")
                        else:
                            prList.append("mtDNA : 1")
                    else:
                        prList.append("mtDNA : N/A")
                # Other inheritance patterns (autosomal dom./rec.,sex-linked (X) dom./rec., sporadic) connected with nuclear mtDNA disorders are not considered here. See the relevant sections: autosomal dom./rec. and sex-linked (X) dom./rec.
                # hz == "F hom" or hz == "M hom" or hz == "F het" or hz == "M het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het":
                else:
                    prList.append("mtDNA : N/A")
            elif pd[i].find("Y-chromosomal") != -1:
                # Paternal inheritance patttern.
                if hz == "":
                    if father == 1 and mother == 0:
                        if sex == "boy":
                            if status == 1:
                                prList.append("Y-chromosomal : 1.0")
                            else:
                                prList.append("Y-chromosomal : 0")
                        elif sex == "girl":
                            if status == 1:
                                prList.append("Y-chromosomal : 0")
                            else:
                                prList.append("Y-chromosomal : 1.0")
                    elif father == 0 and mother == 0:
                        if sex == "boy":
                            if status == 1:
                                prList.append("Y-chromosomal : 0")
                            else:
                                prList.append("Y-chromosomal : 1")
                        elif sex == "girl":
                            if status == 1:
                                prList.append("Y-chromosomal : 0")
                            else:
                                prList.append("Y-chromosomal : 1.0")
                    else:
                        prList.append("Y-chromosomal : N/A")
                    # hz == "F hom" or hz == "M hom" or hz == "F het" or hz == "M het" or hz == "M&F hom" or hz == "F&M hom" or "M&F het" or "F&M het":
                else:
                    prList.append("Y-chromosomal : N/A")

        if status == 1:
            st = "affected"
        else:
            st = "healthy"
        print("Probability of a " + str(sex) + " being " + st + ":")
        for i in range(len(prList)):
            print(prList[i])

    def tbl(self, dr=""):
        pd = self.pdsolver(dr)
        if pd == "Invalid input.":
            return "Invalid input."
        for i in range(1, 27):
            if pd == "? (" + str(i) + ")":
                return "?"
        if pd == "? (sl)" or pd == "? (ad)":
            return "?"
        else:
            for i in range(len(pd)):
                print(pd[i])


# Created by J
