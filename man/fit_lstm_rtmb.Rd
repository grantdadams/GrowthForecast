% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_lstm_rtmb.r
\name{fit_lstm_rtmb}
\alias{fit_lstm_rtmb}
\title{Function to fit LSTM recurrent neural net in RTMB}
\usage{
fit_lstm_rtmb(
  data,
  nhidden_layer = 3,
  hidden_dim = 5,
  input_par = NULL,
  n_proj_years = 2,
  last_year = NA
)
}
\arguments{
\item{data}{data.frame with value, year, and age}

\item{nhidden_layer}{number of hidden layers (this is on a per obs basis)}

\item{hidden_dim}{dimension of hidden layers (nhidden_layer + 2 for in/out)}
}
\description{
Function to fit LSTM recurrent neural net in RTMB
}
