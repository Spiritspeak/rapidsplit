#' Approach-Avoidance Task examining approach bias to different foods
#'
#' This data originates from an approach-avoidance task examining approach bias towards food.
#' Participants responded to the stimulus category (food or object) by pulling or pushing a joystick.
#' Instructions were flipped from one block to the next.
#'
#' @docType data
#'
#' @usage data(foodAAT)
#'
#' @format An object of class \code{"data.frame"}
#' 
#' @details
#' * subjectid: Participant ID
#' * stimid: Stimulus ID
#' * is_pull: Whether the trial required an approach response (1) or an avoid response (0)
#' * is_target: Whether the trial featured a food stimulus (1) or an object stimulus (0)
#' * error: Whether the response was incorrect (1) or correct (0)
#' * RT: The response initiation time
#' * FullRT: The time from stimulus onset to response completion
#' * trialnum: The trial number
#' * blocknum: The block number
#' * palatability: The participant's palatability rating for the stimulus (foods only)
#' * valence: The participant's valence rating for the stimulus
#' * FCQS_2_craving: The participant's FCQS state food craving score at time of testing
#' * FCQS_2_hunger: The participant's FCQS state hunger score at time of testing
#'
#' @keywords datasets
#'
#' @references Lender, A., Meule, A., Rinck, M., Brockmeyer, T., & Blechert, J. (2018). 
#' Measurement of food-related approach–avoidance biases: 
#' Larger biases when food stimuli are task relevant. Appetite, 125, 42–47. 
#' \href{https://doi.org/10.1016/j.appet.2018.01.032}{10.1016/j.appet.2018.01.032}
#' 
#' @source \href{https://doi.org/10.1016/j.appet.2018.01.032}{Original study}
#'
"foodAAT"

# library(magrittr)
# library(dplyr)
# load("./../alldata.rda")
# foodAAT<-datasets$relfoo
# foodAAT %<>% select(subjectid,stimid,is_pull,is_target,error,RT,FullRT,trialnum,blocknum,
#                     palatability,valence,FCQS_2_craving,FCQS_2_hunger)
# save(foodAAT,file="./data/foodAAT.RData")
