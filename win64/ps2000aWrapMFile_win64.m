function [methodinfo,structs,enuminfo,ThunkLibName]=ps2000aWrapMFile_win64
%PS2000AWRAPMFILE_WIN64 Create structures to define interfaces found in 'ps2000aWrap'.

%This function was generated by loadlibrary.m parser version  on Wed Jun 27 14:24:03 2018
%perl options:'ps2000aWrap.i -outfile=ps2000aWrapMFile_win64.m -thunkfile=ps2000aWrap_thunk_pcwin64.c -header=ps2000aWrap.h'
ival={cell(1,0)}; % change 0 to the actual number of functions to preallocate the data.
structs=[];enuminfo=[];fcnNum=1;
fcns=struct('name',ival,'calltype',ival,'LHS',ival,'RHS',ival,'alias',ival,'thunkname', ival);
MfilePath=fileparts(mfilename('fullpath'));
ThunkLibName=fullfile(MfilePath,'ps2000aWrap_thunk_pcwin64');
% extern PICO_STATUS __stdcall RunBlock ( int16_t handle , int32_t preTriggerSamples , int32_t postTriggerSamples , uint32_t timebase , uint32_t segmentIndex ); 
fcns.thunkname{fcnNum}='uint32int16int32int32uint32uint32Thunk';fcns.name{fcnNum}='RunBlock'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'int32', 'int32', 'uint32', 'uint32'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall GetStreamingLatestValues ( int16_t handle ); 
fcns.thunkname{fcnNum}='uint32int16Thunk';fcns.name{fcnNum}='GetStreamingLatestValues'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16'};fcnNum=fcnNum+1;
% extern int32_t __stdcall AvailableData ( int16_t handle , uint32_t * startIndex ); 
fcns.thunkname{fcnNum}='int32int16voidPtrThunk';fcns.name{fcnNum}='AvailableData'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int16', 'uint32Ptr'};fcnNum=fcnNum+1;
% extern int16_t __stdcall AutoStopped ( int16_t handle ); 
fcns.thunkname{fcnNum}='int16int16Thunk';fcns.name{fcnNum}='AutoStopped'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16'};fcnNum=fcnNum+1;
% extern int16_t __stdcall IsReady ( int16_t handle ); 
fcns.thunkname{fcnNum}='int16int16Thunk';fcns.name{fcnNum}='IsReady'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16'};fcnNum=fcnNum+1;
% extern int16_t __stdcall IsTriggerReady ( int16_t handle , uint32_t * triggeredAt ); 
fcns.thunkname{fcnNum}='int16int16voidPtrThunk';fcns.name{fcnNum}='IsTriggerReady'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16', 'uint32Ptr'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall ClearTriggerReady ( int16_t handle ); 
fcns.thunkname{fcnNum}='uint32int16Thunk';fcns.name{fcnNum}='ClearTriggerReady'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall SetTriggerConditions ( int16_t handle , uint32_t * conditionsArray , int16_t nConditions ); 
fcns.thunkname{fcnNum}='uint32int16voidPtrint16Thunk';fcns.name{fcnNum}='SetTriggerConditions'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'uint32Ptr', 'int16'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall SetTriggerProperties ( int16_t handle , int32_t * propertiesArray , int16_t nProperties , int32_t autoTrig ); 
fcns.thunkname{fcnNum}='uint32int16voidPtrint16int32Thunk';fcns.name{fcnNum}='SetTriggerProperties'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'int32Ptr', 'int16', 'int32'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall SetPulseWidthQualifier ( int16_t handle , uint32_t * pwqConditionsArray , int16_t nConditions , uint32_t direction , uint32_t lower , uint32_t upper , uint32_t type ); 
fcns.thunkname{fcnNum}='uint32int16voidPtrint16uint32uint32uint32uint32Thunk';fcns.name{fcnNum}='SetPulseWidthQualifier'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'uint32Ptr', 'int16', 'uint32', 'uint32', 'uint32', 'uint32'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall setChannelCount ( int16_t handle ); 
fcns.thunkname{fcnNum}='uint32int16Thunk';fcns.name{fcnNum}='setChannelCount'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall setEnabledChannels ( int16_t handle , int16_t * enabledChannels ); 
fcns.thunkname{fcnNum}='uint32int16voidPtrThunk';fcns.name{fcnNum}='setEnabledChannels'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'int16Ptr'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall setEnabledDigitalPorts ( int16_t handle , int16_t * enabledDigitalPorts ); 
fcns.thunkname{fcnNum}='uint32int16voidPtrThunk';fcns.name{fcnNum}='setEnabledDigitalPorts'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'int16Ptr'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall setAppAndDriverBuffers ( int16_t handle , int16_t channel , int16_t * appBuffer , int16_t * driverBuffer , int32_t bufferLength ); 
fcns.thunkname{fcnNum}='uint32int16int16voidPtrvoidPtrint32Thunk';fcns.name{fcnNum}='setAppAndDriverBuffers'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'int16', 'int16Ptr', 'int16Ptr', 'int32'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall setMaxMinAppAndDriverBuffers ( int16_t handle , int16_t channel , int16_t * appMaxBuffer , int16_t * appMinBuffer , int16_t * driverMaxBuffer , int16_t * driverMinBuffer , int32_t bufferLength ); 
fcns.thunkname{fcnNum}='uint32int16int16voidPtrvoidPtrvoidPtrvoidPtrint32Thunk';fcns.name{fcnNum}='setMaxMinAppAndDriverBuffers'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'int16', 'int16Ptr', 'int16Ptr', 'int16Ptr', 'int16Ptr', 'int32'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall setAppAndDriverDigiBuffers ( int16_t handle , int16_t digiPort , int16_t * appDigiBuffer , int16_t * driverDigiBuffer , int32_t bufferLength ); 
fcns.thunkname{fcnNum}='uint32int16int16voidPtrvoidPtrint32Thunk';fcns.name{fcnNum}='setAppAndDriverDigiBuffers'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'int16', 'int16Ptr', 'int16Ptr', 'int32'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall setMaxMinAppAndDriverDigiBuffers ( int16_t handle , int16_t digiPort , int16_t * appMaxDigiBuffer , int16_t * appMinDigiBuffer , int16_t * driverMaxDigiBuffer , int16_t * driverMinDigiBuffer , int32_t bufferLength ); 
fcns.thunkname{fcnNum}='uint32int16int16voidPtrvoidPtrvoidPtrvoidPtrint32Thunk';fcns.name{fcnNum}='setMaxMinAppAndDriverDigiBuffers'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='uint32'; fcns.RHS{fcnNum}={'int16', 'int16', 'int16Ptr', 'int16Ptr', 'int16Ptr', 'int16Ptr', 'int32'};fcnNum=fcnNum+1;
structs.tWrapBufferInfo.members=struct('driverBuffers', 'int16PtrPtr', 'appBuffers', 'int16PtrPtr', 'bufferLengths', 'int32#4', 'driverDigiBuffers', 'int16PtrPtr', 'appDigiBuffers', 'int16PtrPtr', 'digiBufferLengths', 'int32#2');
enuminfo.enPS2000AWrapDigitalPortIndex=struct('PS2000A_WRAP_DIGITAL_PORT0',0,'PS2000A_WRAP_DIGITAL_PORT1',1);
methodinfo=fcns;