import Forstat

#Kit name
# readdatafile = Forstat.read_data_excel(kit_name="GlobalFiler",
#                                        file_path="C:\\Users\\peter\\Downloads\\CG_final_data__1_ (1).xlsx",
#                                        consider_multialleles=False
#                                        )
#
readdatafile = Forstat.read_genepop_shorttandemrepeats(kit_name='Goat_TestKit',
                                                       file_path='C:\\Users\\peter\\Downloads\\GOATexample (4).gen',
                                                       digit_format = 3,
                                                       consider_multialleles=False
                                                       )

#If multiallele is True it will determine metrics for indibviduals with multiallelic loci, if False only biallelic loci will be analysed
# #TODO 2. add a user input to ask what analysis should be done analysis will determine the marker type most likely.
#So far we only have autosomalSTRs so we need more applications before this can happen

Forstat.autosomalSTR_population_forensic_metrics(readdatafile)
