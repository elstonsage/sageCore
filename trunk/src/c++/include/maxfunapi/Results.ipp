namespace SAGE   {
namespace MAXFUN {

inline const ParameterMgr & Results::getParameterMgr() const { return my_ParameterMgr; }

inline int Results::setSequenceName      (const std::string & sequence_name)    { my_sequence_name = sequence_name;    return 0; }
inline int Results::setFunctionName      (const std::string & function_name)    { my_function_name = function_name;    return 0; }

inline int Results::setValidInitialValue (bool   valid)   { my_ValidInitialValue = valid;   return 0; }
inline int Results::setWasSkipped        (bool   skipped) { my_WasSkipped        = skipped; return 0; }

inline const std::string & Results::getSequenceName                           () const { return my_sequence_name;     }
inline const std::string & Results::getFunctionName                           () const { return my_function_name;     }

inline double Results::getFinalFunctionValue             () const { return my_FinalFunctionValue;              }
inline int    Results::getExitFlag                       () const { return my_ExitFlag;                        }
inline double Results::getMaximumDifference              () const { return my_MaximumDifference;               }
inline double Results::getERM                            () const { return my_ERM;                             }
inline double Results::getChangeInFunctionValue          () const { return my_ChangeInFunctionValue;           }
inline double Results::getPreviousFunctionValue          () const { return my_PreviousFunctionValue;           }
inline double Results::getGradientNorm                   () const { return my_GradientNorm;                    }
inline double Results::getNewtonGradientNorm             () const { return my_NewtonGradientNorm;              }
inline int    Results::getFirstDerivApproxIndic          () const { return my_FirstDerivApproxIndic;           }
inline int    Results::getAgeOfGradientVector            () const { return my_AgeOfGradientVector;             }
inline int    Results::getGradientStatusFlag             () const { return my_GradientStatusFlag;              }
inline int    Results::getTerminationAtConstraint        () const { return my_TerminationAtConstraint;         }
inline int    Results::getIterations                     () const { return my_Iterations;                      }
inline int    Results::getAgeOfCovMatrix                 () const { return my_AgeOfCovMatrix;                  }
inline int    Results::getCovMatrixStatus                () const { return my_CovMatrixStatus;                 }
inline int    Results::getNumOfConvergedIndependentParams() const { return my_NumOfBoundConvergedIndependentParams; }
inline int    Results::getNumOfDependentParams           () const { return my_NumOfDependentParams;            }
inline int    Results::getNumOfEstimatedParams           () const { return my_NumOfEstimatedParams;            }
inline int    Results::getNumOfIndependentParams         () const { return my_NumOfIndependentParams;          }
inline int    Results::getNumOfVaryingIndependentParams  () const { return my_NumOfVaryingIndependentParams;   }
inline int    Results::getNumOfTrialSearch               () const { return my_NumOfTrialSearch;                }
inline double Results::getStepSizeChange                 () const { return my_StepSizeChange;                  }
inline bool   Results::getValidInitialValue              () const { return my_ValidInitialValue;               }
inline bool   Results::getWasSkipped                     () const { return my_WasSkipped;                      }
                    
inline int    Results::getAugvStatus                     () const { return my_AugvStatus;                      }
} // end namespace MAXFUN
} // end namespace SAGE
