import sys

# arg 1: caller_result.vcf (reslt of the SV caller)
# arg 2: sv_typer_result.vcf

def main():
    caller_genotypes = {}

    # get all SV caller results
    with open(sys.argv[1]) as caller_result_list:
        for line in caller_result_list:
            if line.startswith("#"):
                continue
            tab_split = line.strip().split("\t")
            chr = tab_split[0]
            position = tab_split[1]
            id = chr+"_"+position
            gt = tab_split[9].split(":")[0]

            if gt != "./.":
                caller_genotypes[id] = gt
                
    #print(caller_genotypes)

    # get all SVtyper results
    with open(sys.argv[2]) as svtyper_result_list:
        for line in svtyper_result_list:
            if line.startswith("#"):
                continue

            tab_split = line.strip().split("\t")
            chr = tab_split[0]
            position = tab_split[1]
            id = chr+"_"+position
            sample = tab_split[9].split(":")
            svtyper_gt = sample[0]
            

            if (not id in caller_genotypes):
                print "\t".join(tab_split)
                continue

            if svtyper_gt == "./.":
                sample[0] = caller_genotypes[id]
                tab_split[9] = ":".join(sample)
                print "\t".join(tab_split)
                continue

            # user the SVtyper GT, if it's available
            print "\t".join(tab_split)

            # "Avergage" the GT calls
            # gts = [svtyper_gt, caller_genotypes[id]]

            # het = 0
            # hom = 0
            # ref = 0
            # # counts support for het/hom/ref
            # for gt in gts:
            #     if "0/1" == gt or "1/1" == gt or "./1" == gt or "0/0" == gt:
            #         if "0/1" == gt or "./1" == gt:
            #             het += 1
            #         if "1/1" == gt:
            #             hom += 1
            #         if "0/0" == gt:
            #             ref += 1

            # if het == 0 and hom == 0:
            #     if ref > 0:
            #         sample[0] = "0/0"
            #     else:
            #         sample[0] = "./."
            # elif hom > het:
            #     sample[0] = "1/1"
            # else:
            #     sample[0] = "0/1"

            # tab_split[9] = ":".join(sample)
            # print("\t".join(tab_split))

main()
