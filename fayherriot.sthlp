{smcl}
{* *! version 1.0.0  31December2017}{...}
{cmd:help fayherriot}
{hline}

{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:fayherriot} {hline 1} Fits area level Fay-Herriot model to obtain small area estimates and mean squared error of estimates}{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 23 2}
{opt fayherriot} {varlist} {ifin} {cmd:,}
{opt RE:var(varname)}
[{opt FH:predict(newvarname)}
{opt FHSE:predict(newvarname)}
{opt FHCV:predict(newvarname)}
{opt DSE:predict(newvarname)}
{opt DCV:predict(newvarname)}
{opt AREA:predict(newvarname)}
{opt OUT:sample}]

{title:Description}

{pstd}
{cmd:fayherriot} Supports Fay-Herriot's small area estimation methods, with EBLUP 

{title:Options}

{phang}
{opt RE:var(varname)} Variable containing dependent variable's design-based sampling variance

{phang}
{opt FH:predict(newvarname)} New variable name for predicted Fay Herriot small area estimate

{phang}
{opt FHSE:predict(newvarname)} New variable name for predicted Fay Herriot small area estimate standard error

{phang}
{opt FHCV:predict(newvarname)} New variable name for predicted Fay Herriot small area estimate coefficient of variation

{phang}
{opt DSE:predict(newvarname)} New variable name for predicted Fay Herriot direct estimate standard error 

{phang}
{opt DCV:predict(newvarname)} New variable name for predicted Fay Herriot direct estimate coefficient of variation

{phang}
{opt AREA:predict(newvarname)} New variable name for predicted Fay Herriot area effects

{phang}
{opt OUT:sample} Requests that predictions be made for out of sample observations (requires any of the options fhpredict, fhsepredict, fhcvpredict to be specified )

{phang}
{opt NONEG:ative} Requests that negative predictions be set to 0.

{title:Example}
fayherriot yield hh_f hh_size, revar(vd) fhpredict(fh) fhse(fhse) outsample


{title:Authors}

{pstd}
Paul Corral{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
pcorralrodas@worldbank.org{p_end}

{pstd}
William Seitz{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
wseitz@worldbank.org{p_end}

{title:References}

{pstd}
Chandra, H., Sud, U. C., & Gupta, V. K. (2013). Small Area Estimation under Area Level Model Using R Software.

{pstd}
Fay III, R. E., & Herriot, R. A. (1979). Estimates of income for small places: an application of James-Stein procedures to census data. Journal of the American Statistical Association, 74(366a), 269-277.

