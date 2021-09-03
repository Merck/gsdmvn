#  Copyright (c) 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
#
#  This file is part of the gsdmvn program.
#
#  gsdmvn is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# These global variables are declared to eliminate associated R cmd check warnings.
# There is no other identified functional impact of these global declarations.
# These are column names used in input and/or output tibbles in the package

utils::globalVariables(
  c(
    'Events',
    'Probability',
    'Stratum',
    'Time',
    'AHR',
    'N',
    'Z',
    'dropoutRate',
    'failRate',
    'duration',
    'rate',
    'theta',
    'info',
    'info0',
    'Bound',
    'Analysis',
    'h'
  )
)
