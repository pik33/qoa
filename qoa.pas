unit qoa;

{$mode objfpc}{$H+}

interface

uses
 SysUtils;

// Pascal translation by Piotr Kardasz, pik33@o2.pl

{

Copyright (c) 2023, Dominic Szablewski - https://phoboslab.org
SPDX-License-Identifier: MIT

QOA - The "Quite OK Audio" format for fast, lossy audio compression


-- Data Format

QOA encodes pulse-code modulated (PCM) audio data with up to 255 channels,
sample rates from 1 up to 16777215 hertz and a bit depth of 16 bits.

The compression method employed in QOA is lossy; it discards some information
from the uncompressed PCM data. For many types of audio signals this compression
is "transparent", i.e. the difference from the original file is often not
audible.

QOA encodes 20 samples of 16 bit PCM data into slices of 64 bits. A single
sample therefore requires 3.2 bits of storage space, resulting in a 5x
compression (16 / 3.2).

A QOA file consists of an 8 byte file header, followed by a number of frames.
Each frame contains an 8 byte frame header, the current 16 byte en-/decoder
state per channel and 256 slices per channel. Each slice is 8 bytes wide and
encodes 20 samples of audio data.

All values, including the slices, are big endian. The file layout is as follows:

struct {
	struct {
		char     magic[4];         // magic bytes "qoaf"
		uint32_t samples;          // samples per channel in this file
	} file_header;

	struct {
		struct {
			uint8_t  num_channels; // no. of channels
			uint24_t samplerate;   // samplerate in hz
			uint16_t fsamples;     // samples per channel in this frame
			uint16_t fsize;        // frame size (includes this header)
		} frame_header;

		struct {
			int16_t history[4];    // most recent last
			int16_t weights[4];    // most recent last
		} lms_state[num_channels];

		qoa_slice_t slices[256][num_channels];

	} frames[ceil(samples / (256 * 20))];
} qoa_file_t;

Each `qoa_slice_t` contains a quantized scalefactor `sf_quant` and 20 quantized
residuals `qrNN`:

.- QOA_SLICE -- 64 bits, 20 samples --------------------------/  /------------.
|        Byte[0]         |        Byte[1]         |  Byte[2]  \  \  Byte[7]   |
| 7  6  5  4  3  2  1  0 | 7  6  5  4  3  2  1  0 | 7  6  5   /  /    2  1  0 |
|------------+--------+--------+--------+---------+---------+-\  \--+---------|
|  sf_quant  |  qr00  |  qr01  |  qr02  |  qr03   |  qr04   | /  /  |  qr19   |
`-------------------------------------------------------------\  \------------`

Each frame except the last must contain exactly 256 slices per channel. The last
frame may contain between 1 .. 256 (inclusive) slices per channel. The last
slice (for each channel) in the last frame may contain less than 20 samples; the
slice still must be 8 bytes wide, with the unused samples zeroed out.

Channels are interleaved per slice. E.g. for 2 channel stereo:
slice[0] = L, slice[1] = R, slice[2] = L, slice[3] = R ...

A valid QOA file or stream must have at least one frame. Each frame must contain
at least one channel and one sample with a samplerate between 1 .. 16777215
(inclusive).

If the total number of samples is not known by the encoder, the samples in the
file header may be set to 0x00000000 to indicate that the encoder is
"streaming". In a streaming context, the samplerate and number of channels may
differ from frame to frame. For static files (those with samples set to a
non-zero value), each frame must have the same number of channels and same
samplerate.

Note that this implementation of QOA only handles files with a known total
number of samples.

A decoder should support at least 8 channels. The channel layout for channel
counts 1 .. 8 is:

	1. Mono
	2. L, R
	3. L, R, C
	4. FL, FR, B/SL, B/SR
	5. FL, FR, C, B/SL, B/SR
	6. FL, FR, C, LFE, B/SL, B/SR
	7. FL, FR, C, LFE, B, SL, SR
	8. FL, FR, C, LFE, BL, BR, SL, SR

QOA predicts each audio sample based on the previously decoded ones using a
"Sign-Sign Least Mean Squares Filter" (LMS). This prediction plus the
dequantized residual forms the final output sample.

}

const
 QOA_MIN_FILESIZE=16;
 QOA_MAX_CHANNELS=8;

 QOA_SLICE_LEN=20;
 QOA_SLICES_PER_FRAME=256;
 QOA_FRAME_LEN=(QOA_SLICES_PER_FRAME * QOA_SLICE_LEN);
 QOA_LMS_LEN=4;
 QOA_MAGIC=$716f6166; //* 'qoaf' */



type qoa_lms_t=record
	history:array[0..QOA_LMS_LEN-1] of integer;
	weights:array[0..QOA_LMS_LEN-1] of integer;
        end;
type Pqoa_lms_t=^qoa_lms_t;

type qoa_desc=record
	channels:cardinal;
	samplerate:cardinal;
	samples:cardinal;
        lms:array[0..QOA_MAX_CHANNELS-1] of qoa_lms_t;
        error:double;
	end;

type Pqoa_desc=^qoa_desc;


function qoa_encode_header(qoa:Pqoa_desc; bytes:Pbyte):cardinal;
function qoa_encode_frame(sample_data:Psmallint;qoa:Pqoa_desc; frame_len:cardinal;bytes:Pbyte):cardinal;
function qoa_encode(sample_data:Psmallint;qoa:Pqoa_desc; out_len:Pcardinal):pointer;

function qoa_max_frame_size(qoa:Pqoa_desc):cardinal;
function qoa_decode_header(bytes:Pbyte; size:integer; qoa:Pqoa_desc):cardinal;
function qoa_decode_frame(bytes:Pbyte; size:cardinal; qoa:Pqoa_desc; sample_data:Psmallint; frame_len:Pcardinal):cardinal;
function qoa_decode(bytes:Pbyte; size: integer; qoa:Pqoa_desc):Psmallint;
function qoa_write(filename:Pchar; sample_data:Psmallint; qoa:Pqoa_desc):integer;
function qoa_read(filename:Pchar; qoa:Pqoa_desc):pointer;

implementation



function QOA_MALLOC(sz:cardinal):pointer; inline;
begin result:=getmem(sz); end;
procedure QOA_FREE(p:pointer); inline;
begin freemem(p);end;

function QOA_FRAME_SIZE(channels, slices:integer):integer;

begin result:=(8 + QOA_LMS_LEN * 4 * channels + 8 * slices * channels); end;


type qoa_uint64_t=uint64;


{The quant_tab provides an index into the dequant_tab for residuals in the
range of -8 .. 8. It maps this range to just 3bits and becomes less accurate at
the higher end. Note that the residual zero is identical to the lowest positive
value. This is mostly fine, since the qoa_div() function always rounds away
from zero.}

type Tquoa_quant_tab=array[0..16] of integer;
const qoa_quant_tab:Tquoa_quant_tab =
(	7, 7, 7, 5, 5, 3, 3, 1, ///* -8..-1 */
	0,                      ///*  0     */
	0, 2, 2, 4, 4, 6, 6, 6  ///*  1.. 8 */
);


{We have 16 different scalefactors. Like the quantized residuals these become
less accurate at the higher end. In theory, the highest scalefactor that we
would need to encode the highest 16bit residual is (2**16)/8 = 8192. However we
rely on the LMS filter to predict samples accurately enough that a maximum
residual of one quarter of the 16 bit range is sufficient. I.e. with the
scalefactor 2048 times the quant range of 8 we can encode residuals up to 2**14.

The scalefactor values are computed as:
scalefactor_tab[s] <- round(pow(s + 1, 2.75)) }

type Tqoa_tab=array[0..15] of integer;
const qoa_scalefactor_tab:Tqoa_tab = (
	1, 7, 21, 45, 84, 138, 211, 304, 421, 562, 731, 928, 1157, 1419, 1715, 2048
);


{The reciprocal_tab maps each of the 16 scalefactors to their rounded
reciprocals 1/scalefactor. This allows us to calculate the scaled residuals in
the encoder with just one multiplication instead of an expensive division. We
do this in .16 fixed point with integers, instead of floats.

The reciprocal_tab is computed as:
reciprocal_tab[s] <- ((1<<16) + scalefactor_tab[s] - 1) / scalefactor_tab[s] }

const qoa_reciprocal_tab:Tqoa_tab = (
	65536, 9363, 3121, 1457, 781, 475, 311, 216, 156, 117, 90, 71, 57, 47, 39, 32
);


{The dequant_tab maps each of the scalefactors and quantized residuals to
their unscaled & dequantized version.

Since qoa_div rounds away from the zero, the smallest entries are mapped to 3/4
instead of 1. The dequant_tab assumes the following dequantized values for each
of the quant_tab indices and is computed as:
float dqt[8] = {0.75, -0.75, 2.5, -2.5, 4.5, -4.5, 7, -7};
dequant_tab[s][q] <- round_ties_away_from_zero(scalefactor_tab[s] * dqt[q])

The rounding employed here is "to nearest, ties away from zero",  i.e. positive
and negative values are treated symmetrically.
}

type Tqoa_dequant_tab=array[0..15,0..7] of integer;
const qoa_dequant_tab:Tqoa_dequant_tab = (
	(   1,    -1,    3,    -3,    5,    -5,     7,     -7),
	(   5,    -5,   18,   -18,   32,   -32,    49,    -49),
	(  16,   -16,   53,   -53,   95,   -95,   147,   -147),
	(  34,   -34,  113,  -113,  203,  -203,   315,   -315),
	(  63,   -63,  210,  -210,  378,  -378,   588,   -588),
	( 104,  -104,  345,  -345,  621,  -621,   966,   -966),
	( 158,  -158,  528,  -528,  950,  -950,  1477,  -1477),
	( 228,  -228,  760,  -760, 1368, -1368,  2128,  -2128),
	( 316,  -316, 1053, -1053, 1895, -1895,  2947,  -2947),
	( 422,  -422, 1405, -1405, 2529, -2529,  3934,  -3934),
	( 548,  -548, 1828, -1828, 3290, -3290,  5117,  -5117),
	( 696,  -696, 2320, -2320, 4176, -4176,  6496,  -6496),
	( 868,  -868, 2893, -2893, 5207, -5207,  8099,  -8099),
	(1064, -1064, 3548, -3548, 6386, -6386,  9933,  -9933),
	(1286, -1286, 4288, -4288, 7718, -7718, 12005, -12005),
	(1536, -1536, 5120, -5120, 9216, -9216, 14336, -14336)
);


{ The Least Mean Squares Filter is the heart of QOA. It predicts the next
sample based on the previous 4 reconstructed samples. It does so by continuously
adjusting 4 weights based on the residual of the previous prediction.

The next sample is predicted as the sum of (weight[i] * history[i]).

The adjustment of the weights is done with a "Sign-Sign-LMS" that adds or
subtracts the residual to each weight, based on the corresponding sample from
the history. This, surprisingly, is sufficient to get worthwhile predictions.

This is all done with fixed point integers. Hence the right-shifts when updating
the weights and calculating the prediction. }

function qoa_lms_predict(lms:Pqoa_lms_t):integer; inline;

var prediction:integer=0;
    i:integer;

begin
for i := 0 to QOA_LMS_LEN do prediction += lms^.weights[i] * lms^.history[i];
result:=prediction shr 13;
end;


procedure qoa_lms_update(lms:Pqoa_lms_t; sample,residual:integer); inline;

var delta,i:integer;

begin
delta := residual shr 4;
for i := 0 to QOA_LMS_LEN-1 do
  if lms^.history[i] < 0 then lms^.weights[i] -= delta else lms^.weights[i] +=delta;
for i := 0 to QOA_LMS_LEN-2 do lms^.history[i] := lms^.history[i+1];
lms^.history[QOA_LMS_LEN-1] := sample;
end;


{qoa_div() implements a rounding division, but avoids rounding to zero for
small numbers. E.g. 0.1 will be rounded to 1. Note that 0 itself still
returns as 0, which is handled in the qoa_quant_tab[].
qoa_div() takes an index into the .16 fixed point qoa_reciprocal_tab as an
argument, so it can do the division with a cheaper integer multiplication. }

function qoa_div(v,scalefactor:integer):integer;inline;

var reciprocal,n,v1,v2,n1,n2:integer;

begin
reciprocal := qoa_reciprocal_tab[scalefactor];
n := (v * reciprocal + (1 shl 15)) shr 16;
if v>0 then v1:=1 else v1:=0;
if v<0 then v2:=1 else v2:=0;
if n>0 then n1:=1 else n1:=0;
if n<0 then n2:=1 else n2:=0;

n := n + (v1-v2) - (n1-n2); //* round away from 0 */
result:=n;
end;

function qoa_clamp(v,min,max:integer):integer; inline;

begin
result:=v;
if (v < min) then result:=min;
if (v > max) then result:=max;
end;

{This specialized clamp function for the signed 16 bit range improves decode
performance quite a bit. The extra if() statement works nicely with the CPUs
branch prediction as this branch is rarely taken.}


function qoa_clamp_s16(v:integer):integer; inline;

begin
        result:=v;
	if cardinal(v + 32768) > 65535 then
                begin
		if (v < -32768) then result:=-32768;
		if (v >  32767) then result:=32767;
                end;
end;


function qoa_read_u64(bytes:Pbyte;p:Pcardinal): qoa_uint64_t; inline;

begin
bytes += p^;
p^ += 8;
result:=(qoa_uint64_t(bytes[0]) shl 56) or (qoa_uint64_t(bytes[1]) shl 48) or
		(qoa_uint64_t(bytes[2]) shl 40) or (qoa_uint64_t(bytes[3]) shl 32) or
		(qoa_uint64_t(bytes[4]) shl 24) or (qoa_uint64_t(bytes[5]) shl 16) or
		(qoa_uint64_t(bytes[6]) shl  8) or (qoa_uint64_t(bytes[7]) shl  0);
end;

procedure qoa_write_u64(v:qoa_uint64_t; bytes:Pbyte; p:Pcardinal); inline;

begin
bytes += p^;
p^ += 8;
bytes[0] := (v shr 56) and $ff;
bytes[1] := (v shr 48) and  $ff;
bytes[2] := (v shr 40) and  $ff;
bytes[3] := (v shr 32) and  $ff;
bytes[4] := (v shr 24) and  $ff;
bytes[5] := (v shr 16) and  $ff;
bytes[6] := (v shr 8) and  $ff;
bytes[7] := (v shr 0) and  $ff;
end;


///* -----------------------------------------------------------------------------
//	Encoder */

function qoa_encode_header(qoa:Pqoa_desc; bytes:Pbyte):cardinal;

var p:cardinal;

begin
p := 0;
qoa_write_u64((qoa_uint64_t(QOA_MAGIC) shl 32) or qoa^.samples,bytes, @p);
result:=p;
end;


function qoa_encode_frame(sample_data:Psmallint; qoa:Pqoa_desc; frame_len:cardinal; bytes:Pbyte):cardinal;

var i, c, channels, p, slices,frame_size, sample_index :cardinal;
    slice_len,slice_start,slice_end,best_scalefactor,si2,sfi,scalefactor,si3:integer;
    prev_scalefactor:array[0..QOA_MAX_CHANNELS-1] of integer;
    best_rank,best_error,weights,history,best_slice,slice,current_rank,current_error,error_sq:qoa_uint64_t;
    lms, best_lms:qoa_lms_t;
    sample,predicted,residual,scaled,clamped,quantized,dequantized,reconstructed,si,weights_penalty:integer;
    error:int64;

begin
channels := qoa^.channels;
p := 0;
slices := (frame_len + QOA_SLICE_LEN - 1) div QOA_SLICE_LEN;
frame_size := QOA_FRAME_SIZE(channels, slices);
for i:=0 to QOA_MAX_CHANNELS-1 do prev_scalefactor[i]:=0;

// Write the frame header

qoa_write_u64((
		qoa_uint64_t(qoa^.channels)   shl 56 or
		qoa_uint64_t(qoa^.samplerate) shl 32 or
		qoa_uint64_t(frame_len)       shl 16 or
		qoa_uint64_t(frame_size)
	), bytes, @p);


for c := 0 to channels-1 do
  begin

  // Write the current LMS state */
  weights := 0;
  history := 0;
  for i := 0 to QOA_LMS_LEN-1 do
    begin
    history := (history shl 16) or (qoa^.lms[c].history[i] and $ffff);
    weights := (weights shl 16) or (qoa^.lms[c].weights[i] and $ffff);
    end;
  qoa_write_u64(history, bytes, @p);
  qoa_write_u64(weights, bytes, @p);
  end;

// We encode all samples with the channels interleaved on a slice level.
// E.g. for stereo: (ch-0, slice 0), (ch 1, slice 0), (ch 0, slice 1), ...

si2:=-QOA_SLICE_LEN; // helper var for sample_index for c translated loop
for sample_index := 0 to (frame_len div QOA_SLICE_LEN)-1 do
  begin
  si2+=QOA_SLICE_LEN;
  for c := 0 to channels-1 do
    begin
    slice_len := qoa_clamp(QOA_SLICE_LEN, 0, frame_len - si2);
    slice_start := si2 * channels + c;
    slice_end := (si2 + slice_len) * channels + c;

    // Brute force search for the best scalefactor. Just go through all
    // 16 scalefactors, encode all samples for the current slice and
    // meassure the total squared error.


    best_rank :=  -1; //18446744073709551615;
    best_error := -1; //18446744073709551615;
    best_slice := 0;
    best_scalefactor := 0;
    for sfi := 0 to 15 do
				//* There is a strong correlation between the scalefactors of
				// neighboring slices. As an optimization, start testing
				// the best scalefactor of the previous slice first. */
      begin
      scalefactor := (sfi + prev_scalefactor[c]) mod 16;

				//* We have to reset the LMS state to the last known good one
				//before trying each scalefactor, as each pass updates the LMS
				//state when encoding. */
      lms := qoa^.lms[c];
      slice := scalefactor;
      current_rank := 0;
      current_error := 0;

      si3:=-channels;
      for si := slice_start to (slice_end div channels)-1 do
        begin
        si3+=channels;
        sample := sample_data[si3];
	predicted := qoa_lms_predict(@lms);
        residual := sample - predicted;
	scaled := qoa_div(residual, scalefactor);
	clamped := qoa_clamp(scaled, -8, 8);
	quantized := qoa_quant_tab[clamped + 8];
	dequantized := qoa_dequant_tab[scalefactor][quantized];
	reconstructed := qoa_clamp_s16(predicted + dequantized);


					//* If the weights have grown too large, we introduce a penalty
					//here. This prevents pops/clicks in certain problem cases */

	weights_penalty := ((
	lms.weights[0] * lms.weights[0] +
	lms.weights[1] * lms.weights[1] +
	lms.weights[2] * lms.weights[2] +
	lms.weights[3] * lms.weights[3]
	) shr 18) - $8ff;
	if (weights_penalty < 0) then weights_penalty := 0;
	error := (sample - reconstructed);
	error_sq := error * error;
        current_rank += error_sq + weights_penalty * weights_penalty;
	current_error += error_sq;
	if (current_rank > best_rank) then break;
	qoa_lms_update(@lms, reconstructed, dequantized);
	slice := (slice shl 3) or quantized;
	end;
      if (current_rank < best_rank) then
        begin
	best_rank := current_rank;
	best_error := current_error;
	best_slice := slice;
	best_lms := lms;
	best_scalefactor := scalefactor;
	end;
      end;
    prev_scalefactor[c] := best_scalefactor;
    qoa^.lms[c] := best_lms;
    qoa^.error += best_error;

                       	//* If this slice was shorter than QOA_SLICE_LEN, we have to left-
			//shift all encoded data, to ensure the rightmost bits are the empty
			//ones. This should only happen in the last frame of a file as all
			//slices are completely filled otherwise. */
     best_slice :=best_slice shl ((QOA_SLICE_LEN - slice_len) * 3);
     qoa_write_u64(best_slice, bytes, @p);
    end;
  end;
result:=p;
end;


function qoa_encode(sample_data:Psmallint; qoa:Pqoa_desc; out_len:Pcardinal):pointer;

var  c,p,i,num_frames, num_slices, encoded_size,sample_index,  frame_size:cardinal;
     bytes:Pbyte;
     frame_len:integer;
     frame_samples:Psmallint;

begin
if (qoa^.samples = 0) or (qoa^.samplerate = 0) or (qoa^.samplerate > $ffffff) or
   (qoa^.channels = 0)  or (qoa^.channels > QOA_MAX_CHANNELS)
   then exit(nil);

//* Calculate the encoded size and allocate */

num_frames := (qoa^.samples + QOA_FRAME_LEN-1) div QOA_FRAME_LEN;
num_slices := (qoa^.samples + QOA_SLICE_LEN-1) div QOA_SLICE_LEN;
encoded_size := 8 +                                            { 8 byte file header }
		num_frames * 8 +                               { 8 byte frame headers }
		num_frames * QOA_LMS_LEN * 4 * qoa^.channels + { 4 * 4 bytes lms state per channel }
		num_slices * 8 * qoa^.channels;                { 8 byte slices }

bytes := QOA_MALLOC(encoded_size);
for c := 0 to  qoa^.channels-1 do
  begin
                //* Set the initial LMS weights to {0, 0, -1, 2}. This helps with the
		//prediction of the first few ms of a file. */
  qoa^.lms[c].weights[0] := 0;
  qoa^.lms[c].weights[1] := 0;
  qoa^.lms[c].weights[2] := -(1<<13);
  qoa^.lms[c].weights[3] :=  (1<<14);

		///* Explicitly set the history samples to 0, as we might have some
		//garbage in there. */
  for i := 0 to QOA_LMS_LEN-1 do qoa^.lms[c].history[i] := 0;
  end;

        //* Encode the header and go through all frames */

p := qoa_encode_header(qoa, bytes);
qoa^.error := 0;
frame_len := QOA_FRAME_LEN;

sample_index := 0;
repeat
  frame_len := qoa_clamp(QOA_FRAME_LEN, 0, qoa^.samples - sample_index);
  frame_samples := sample_data + sample_index * qoa^.channels;
  frame_size := qoa_encode_frame(frame_samples, qoa, frame_len, bytes + p);
  p += frame_size;
  sample_index+=frame_len
until sample_index>=qoa^.samples;
out_len^ := p;
result:= bytes;
end;



//  -----------------------------------------------------------------------------
//	Decoder

function qoa_max_frame_size(qoa:Pqoa_desc):cardinal;

begin
result:= QOA_FRAME_SIZE(qoa^.channels, QOA_SLICES_PER_FRAME);
end;

function qoa_decode_header(bytes:Pbyte; size:integer; qoa:PQoa_desc):cardinal;

var p:cardinal;
    file_header,frame_header:qoa_uint64_t;

begin
p := 0;
if (size < QOA_MIN_FILESIZE) then exit(0);
	// Read the file header, verify the magic number ('qoaf') and read the
	// total number of samples. */
file_header := qoa_read_u64(bytes, @p);
if ((file_header >> 32) <> QOA_MAGIC) then exit(0);
qoa^.samples := file_header and $ffffffff;
if (qoa^.samples=0) then exit(0);
	// Peek into the first frame header to get the number of channels and
	//the samplerate. */
frame_header := qoa_read_u64(bytes, @p);
qoa^.channels := (frame_header shr 56) and $0000ff;
qoa^.samplerate := (frame_header shr 32) and $ffffff;
if (qoa^.channels = 0) or (qoa^.samples = 0) or(qoa^.samplerate = 0) then exit(0);
result:= 8;
end;


function qoa_decode_frame(bytes:Pbyte; size:cardinal; qoa:PQoa_desc; sample_data:Psmallint; frame_len:Pcardinal): cardinal;

var p, c, i, channels,samplerate,samples,frame_size,data_size,num_slices,max_total_samples,sample_index,si2:cardinal;
    frame_header, history, weights, slice :qoa_uint64_t;
    scalefactor, slice_start, slice_end, si, predicted, quantized, dequantized, reconstructed:integer;

begin
p := 0;
frame_len^ := 0;
if (size < 8 + QOA_LMS_LEN * 4 * qoa^.channels) then exit(0);
	// Read and verify the frame header */
frame_header := qoa_read_u64(bytes, @p);
channels   := (frame_header shr 56) and $0000ff;
samplerate := (frame_header shr 32) and $ffffff;
samples    := (frame_header shr 16) and $00ffff;
frame_size := (frame_header      ) and $00ffff;

data_size := frame_size - 8 - QOA_LMS_LEN * 4 * channels;
num_slices := data_size div 8;
max_total_samples := num_slices * QOA_SLICE_LEN;
if (channels <> qoa^.channels) or (samplerate <> qoa^.samplerate) or (frame_size > size) or (samples*channels > max_total_samples) then exit(0);
	// Read the LMS state: 4 x 2 bytes history, 4 x 2 bytes weights per channel */
for c := 0 to channels-1 do
  begin
  history := qoa_read_u64(bytes, @p);
  weights := qoa_read_u64(bytes, @p);
  for i := 0 to QOA_LMS_LEN-1 do
    begin
    qoa^.lms[c].history[i] := smallint(history shr 48);
    history := history shl 16;
    qoa^.lms[c].weights[i] := smallint(weights shr 48);
    weights :=weights shl 16;
    end;
  end;
 	// Decode all slices for all channels in this frame */

for si2 := 0 to (samples div QOA_SLICE_LEN)-1 do
  begin
  sample_index:=si2*QOA_SLICE_LEN;
  for c := 0 to channels-1 do
    begin
    slice := qoa_read_u64(bytes, @p);
    scalefactor := (slice shr 60) and $f;
    slice:= slice shl 4;
    slice_start := sample_index * channels + c;
    slice_end := qoa_clamp(sample_index + QOA_SLICE_LEN, 0, samples) * channels + c;
    si:=slice_start;
    repeat
      predicted := qoa_lms_predict(@(qoa^.lms[c]));
      quantized := (slice shr 61) and $7;
      dequantized := qoa_dequant_tab[scalefactor][quantized];
      reconstructed := qoa_clamp_s16(predicted + dequantized);
      sample_data[si] := reconstructed;
      slice:=slice shl 3;
      qoa_lms_update(@(qoa^.lms[c]), reconstructed, dequantized);
      si += channels;
    until si>=slice_end
   end;
  end;
frame_len^ := samples;
result:=p;
end;

function qoa_decode(bytes:Pbyte; size:integer; qoa:Pqoa_desc):Psmallint;

var p:cardinal;
    total_samples:integer;
    sample_ptr, sample_data:Psmallint;
    sample_index, frame_len, frame_size:cardinal;

begin
p := qoa_decode_header(bytes, size, qoa);
if p=0 then exit(nil);

	//* Calculate the required size of the sample buffer and allocate */
total_samples := qoa^.samples * qoa^.channels;
sample_data := QOA_MALLOC(total_samples * sizeof(smallint));
sample_index := 0;
	//* Decode all frames */

repeat
  sample_ptr := sample_data + sample_index * qoa^.channels;
  frame_size := qoa_decode_frame(bytes + p, size - p, qoa, sample_ptr, @frame_len);
  p += frame_size;
  sample_index += frame_len;
until (frame_size=0) or (sample_index >= qoa^.samples);
qoa^.samples := sample_index;
result:=sample_data;
end;



// -----------------------------------------------------------------------------
//File read/write convenience functions */

function qoa_write(filename:Pchar; sample_data:Psmallint; qoa:Pqoa_desc):integer;

var size:cardinal;
    encoded:pbyte;
    f:integer;

begin
f:= filecreate(filename);
//f := fileopen(filename, fmOpenWrite);
if f<=0 then exit(0);
encoded := qoa_encode(sample_data, qoa, @size);
if encoded=nil then begin fileclose(f); exit(0); end;
filewrite(f,encoded^,size);
fileclose(f);
QOA_FREE(encoded);
result:=size
end;

function qoa_read(filename:Pchar; qoa:Pqoa_desc):pointer;

var f:integer;
    size,bytes_read:cardinal;
    data:pbyte;
    sample_data:psmallint;

begin
f := fileopen(filename, fmOpenRead);
if f<=0 then exit(nil);
size:=fileseek(f,0,fsFromEnd);
if size <=0 then begin fileclose(f); exit(nil); end;
fileseek(f,0,fsFromBeginning);
data := QOA_MALLOC(size);
if (data=nil) then begin fileclose(f); exit(nil); end;
bytes_read := fileread(f,data^, size);
fileclose(f);
sample_data := qoa_decode(data, bytes_read, qoa);
QOA_FREE(data);
result:= sample_data;
end;

end.

