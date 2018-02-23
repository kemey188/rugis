// Copyright 2010 Google
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef BASE_COMMANDLINEFLAGS_H
#define BASE_COMMANDLINEFLAGS_H

// #include "gflags/gflags.h"
#include <string>
#include <vector>

// We care a lot about number of bits things take up.  Unfortunately,
// systems define their bit-specific ints in a lot of different ways.
// We use our own way, and have a typedef to get there.
// Note: these commands below may look like "#if 1" or "#if 0", but
// that's because they were constructed that way at ./configure time.
// Look at gflags.h.in to see how they're calculated (based on your config).
#if 1
#include <stdint.h>             // the normal place uint16_t is defined
#endif
#if 1
#include <sys/types.h>          // the normal place u_int16_t is defined
#endif
#if 1
#include <inttypes.h>           // a third place for uint16_t or u_int16_t
#endif



#endif  // BASE_COMMANDLINEFLAGS_H
