/**********************************************************************
  Copyright(c) 2011-2024 Intel Corporation All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************/

/* Altered from igzip.c to retain only helpers referenced by inflate. */

#include <stddef.h>
#include <stdint.h>

#include "igzip_checksums.h"
#include "igzip_lib.h"

uint32_t
adler32_base(uint32_t adler32, uint8_t *start, uint64_t length);

uint32_t
isal_adler32_bam1(uint32_t adler32, const unsigned char *start, uint64_t length)
{
        uint32_t a = adler32 & 0xffff;

        a = (a == ADLER_MOD - 1) ? 0 : a + 1;
        adler32 = adler32_base((adler32 & 0xffff0000) | a, (uint8_t *) start, length);
        a = adler32 & 0xffff;
        a = (a == 0) ? ADLER_MOD - 1 : a - 1;
        return (adler32 & 0xffff0000) | a;
}

void
isal_gzip_header_init(struct isal_gzip_header *gz_hdr)
{
        gz_hdr->text = 0;
        gz_hdr->time = 0;
        gz_hdr->xflags = 0;
        gz_hdr->os = 0xff;
        gz_hdr->extra = NULL;
        gz_hdr->extra_buf_len = 0;
        gz_hdr->extra_len = 0;
        gz_hdr->name = NULL;
        gz_hdr->name_buf_len = 0;
        gz_hdr->comment = NULL;
        gz_hdr->comment_buf_len = 0;
        gz_hdr->hcrc = 0;
        gz_hdr->flags = 0;
}

void
isal_zlib_header_init(struct isal_zlib_header *z_hdr)
{
        z_hdr->info = 0;
        z_hdr->level = 0;
        z_hdr->dict_id = 0;
        z_hdr->dict_flag = 0;
}
