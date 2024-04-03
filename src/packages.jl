# =================================================================
# INCLUDED PACKAGES
# =================================================================

using LinearAlgebra         # standard linear algebra routines

if ~Sys.isapple()
    using MKL                   # use Intel math libs
    ENV["MKL_DYNAMIC"] = false
end

BLAS.set_num_threads(1)     # use BLAS single threaded 
using Parameters            # handle parameters more efficiently
using Printf                # nicer console output
using Random                # for drawing random variables
using Distributions         # for drawing random variables
using Statistics            # for all stat functions
using LaTeXStrings          # used for Latex in plots
using FastGaussQuadrature   # quadrature integration 
using StatsBase             # standard statistics 
using NLsolve               # non linear solver 
using CSV                   # import CSV
using DataFrames            # data frame structures
using Plots.PlotMeasures    # allows to use mm for margins, size etc in plots
using DelimitedFiles        # read csv files 
using Interpolations        # spline interpolation
using Dates                 # handle dates
using XLSX                  # excel files
using FileIO, JLD2          # file storage 
using TimerOutputs          # timing tool
using StaticArrays          # use static arrays


using Plots                 # for figures
using CairoMakie
# using GLMakie

# Set latex scheme for Cairo
MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Bold.otf")
))

# using QuantEcon
