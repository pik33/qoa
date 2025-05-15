unit Unit1;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls;

type

  { TForm1 }

  TForm1 = class(TForm)
    Button1: TButton;
    Button2: TButton;
    OpenDialog1: TOpenDialog;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
  private

  public

  end;

var
  Form1: TForm1;

implementation

{$R *.lfm}

uses qoa;

{
Converter translated to freepascal by pik33@02.pl.
Original C code Copyright (c) 2023, Dominic Szablewski - https://phoboslab.org
SPDX-License-Identifier: MIT
}

// -----------------------------------------------------------------------------
//	WAV reader / writer */

function QOACONV_CHUNK_ID(S:pchar):cardinal;
begin
result:= (cardinal(S[3]) shl 24) or (cardinal(S[2]) shl 16) or (cardinal(S[1]) shl 8) or (cardinal(S[0])) ;
end;

procedure qoaconv_fwrite_u32_le(v,fh:integer);

var buf:array[0..3] of byte;
    buf1:cardinal absolute buf;

begin
buf[0] := $ff and (v      );
buf[1] := $ff and (v shr  8);
buf[2] := $ff and (v shr 16);
buf[3] := $ff and (v shr 24);
filewrite(fh,buf,4);
end;

procedure qoaconv_fwrite_u16_le(v:word; fh:integer);

var buf:array[0..1] of byte;
    wrote:integer;

begin

buf[0] := $ff and (v      );
buf[1] := $ff and (v shr  8);
wrote:=filewrite(fh,buf,2);
// if wrote<>0 then throw write error
end;


function qoaconv_fread_u32_le(fh:integer):cardinal;

var buf:array[0..3] of byte;
    rd:integer;

begin
rd := fileread(fh,buf,4);
// if rd=0 then throw read error
result := (buf[3] shl 24) or (buf[2] shl 16) or (buf[1] shl 8) or buf[0];
end;

function  qoaconv_fread_u16_le(fh: integer):word;

var buf:array[0..1] of byte;
    rd:integer;

begin
rd:=fileread(fh,buf,2);
// if rd=0 then throw read error
result:= (buf[1] shl 8) or buf[0];
end;

function qoaconv_wav_write(path:Pchar; sample_data: PSmallint; desc:Pqoa_desc):integer;

var data_size, ds2,samplerate:cardinal;
    channels,bits_per_sample:word;
    fh:integer;
    s:Pchar;

begin
data_size := desc^.samples * desc^.channels * sizeof(smallint);
samplerate := desc^.samplerate;
bits_per_sample := 16;
channels := desc^.channels;

	//* Lifted from https://www.jonolick.com/code.html - public domain
	//Made endian agnostic using qoaconv_fwrite() */

fh := filecreate(path);
//	QOACONV_ASSERT(fh, "Can't open %s for writing", path);
filewrite(fh,'RIFF',4);
ds2:=data_size+44-8;
qoaconv_fwrite_u32_le(ds2,fh);
s:='WAVEfmt '+#16+#0+#0+#0+#1+#0;
filewrite(fh, s^, 14);
qoaconv_fwrite_u16_le(channels, fh);
qoaconv_fwrite_u32_le(samplerate, fh);
qoaconv_fwrite_u32_le(channels * samplerate * bits_per_sample div 8, fh);
qoaconv_fwrite_u16_le(channels * bits_per_sample div 8, fh);
qoaconv_fwrite_u16_le(bits_per_sample, fh);
filewrite(fh,'data', 4);
qoaconv_fwrite_u32_le(data_size, fh);
filewrite(fh,sample_data^, data_size);
fileclose(fh);
result:= data_size  + 44 - 8;
end;

function qoaconv_wav_read(path:PChar; desc: Pqoa_desc):Psmallint;

var fh:integer;
    buf:array[0..43] of byte;
    buf2: array[0..21] of word absolute buf;
    buf3:array[0..10] of cardinal absolute buf;
    spl:Psmallint;

begin
fh := fileopen(path, fmopenread);
fileread(fh,buf,44);


desc^.samplerate := buf3[6];
desc^.samples := buf3[10] div (buf2[11] * buf2[17] div 8);
desc^.channels := buf2[11];
spl:=getmem(buf3[10]);
fileread(fh,spl^,buf3[10]);
result:=spl;
fileclose(fh);
end;


procedure encode (filename:string);


var desc:qoa_desc;
    sample_data:Psmallint;


begin
sample_data := qoaconv_wav_read(pchar(filename), @desc);
qoa_write(pchar(filename+'.qoa'), sample_data, @desc);
end;


procedure decode  (filename:string);

var desc:qoa_desc;
sample_data:Psmallint;
bytes_written:cardinal;

begin
sample_data := qoa_read(pchar(filename), @desc);
bytes_written := qoaconv_wav_write(pchar(filename+'.wav'), sample_data, @desc);
end;

{ TForm1 }

procedure TForm1.Button1Click(Sender: TObject);
begin
  if opendialog1.execute then encode(opendialog1.filename);
end;

procedure TForm1.Button2Click(Sender: TObject);
begin
  if opendialog1.execute then decode(opendialog1.filename);
end;






end.

