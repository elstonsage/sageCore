#########################################################################
# This data file was generated using the config file 'wedig.cfg'
#########################################################################

pedigree
{
  delimiter_mode = multiple
  delimiters="	"
  individual_missing_value="0"
  sex_code, male="0", female="1", unknown="?"

  pedigree_id   = FAM
  individual_id = ID
  parent_id     = MOTHER
  parent_id     = FATHER
  sex_field     = SEX

  marker = MAR1, delimiter = /, name = MAR1, missing = MISS
  trait = MAR1pheno, missing = MISS
  trait = cov, missing = MISS
  trait = t1, missing = MISS
  trait = t2, missing = MISS
  trait = t3, missing = MISS
}

segreg
{
  trait=t1
  class=D

  transformation
  {
    option=none
  }

  type_mean
  {
    option=one
  }

  mean_cov
  {
    covariate=cov
  }
  resid
  {
    option=arb
    fm=0,fixed=false
    ss=0.1
  }

}
