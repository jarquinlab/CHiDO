# CHiDO is a no-code platform to integrate multi-omics data to build, train and test
# linear mixed models for identifying candidates for desired GxE interactions.
#
# Copyright (C) 2024 Francisco Gonzalez, Diego Jarquin, and Julian Garcia
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Created by: Francisco Gonzalez
# Last updated by: Francisco Gonzalez
# Last updated: 06/09/2024

source("code/functions/utils.R")

### UI instructions ----
### Supporting objects ----

about_md <- markdown::renderMarkdown(file="README.md")

### Main tab ----

about_tab <- tabItem(
  tabName = "about",
  HTML(about_md)
)