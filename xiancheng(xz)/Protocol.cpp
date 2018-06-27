/****************************************Copyright(c)*****************************************************
**                            Shenzhen Yuejiang Technology Co., LTD.
**
**                                 http://www.dobot.cc
**
**--------------File Info---------------------------------------------------------------------------------
** File name:           Protocol.cpp
** Latest modified Date:2016-06-01
** Latest Version:      V1.0.0
** Descriptions:        Protocol interface
**
**--------------------------------------------------------------------------------------------------------
** Created by:          Liu Zhufu
** Created date:        2016-03-14
** Version:             V1.0.0
** Descriptions:
**--------------------------------------------------------------------------------------------------------
*********************************************************************************************************/

#include "Protocol.h"
#include <stdio.h>
#include <string.h>
#include "command.h"



/*********************************************************************************************************
** Protocol buffer definition
*********************************************************************************************************/
#define RAW_TX_BYTE_BUFFER_SIZE    1024
#define RAW_RX_BYTE_BUFFER_SIZE    4096
#define RAW_BYTE_BUFFER_SIZE    4096
#define PACKET_BUFFER_SIZE  8

// UART4
uint8_t gUART4TXRawByteBuffer[RAW_BYTE_BUFFER_SIZE];
uint8_t gUART4RXRawByteBuffer[RAW_BYTE_BUFFER_SIZE];
Packet gUART4TXPacketBuffer[PACKET_BUFFER_SIZE];
Packet gUART4RXPacketBuffer[PACKET_BUFFER_SIZE];

ProtocolHandler gUART4ProtocolHandler;

//add by meng jun
extern int flag_tx;


/*********************************************************************************************************
** Function name:       ProtocolInit
** Descriptions:        Init the protocol buffer etc.
** Input parameters:    None
** Output parameters:   None
** Returned value:      None
*********************************************************************************************************/
void ProtocolInit(void)
{
    // Init UART4 protocol
    RingBufferInit(&gUART4ProtocolHandler.txRawByteQueue, gUART4TXRawByteBuffer, RAW_BYTE_BUFFER_SIZE, sizeof(uint8_t));
    RingBufferInit(&gUART4ProtocolHandler.rxRawByteQueue, gUART4RXRawByteBuffer, RAW_BYTE_BUFFER_SIZE, sizeof(uint8_t));
    RingBufferInit(&gUART4ProtocolHandler.txPacketQueue, gUART4TXPacketBuffer, PACKET_BUFFER_SIZE, sizeof(Packet));
    RingBufferInit(&gUART4ProtocolHandler.rxPacketQueue, gUART4RXPacketBuffer, PACKET_BUFFER_SIZE, sizeof(Packet));

}

/*********************************************************************************************************
** Function name:       ProtocolProcess
** Descriptions:        Process the protocol
** Input parameters:    None
** Output parameters:   None
** Returned value:      None
*********************************************************************************************************/
void ProtocolProcess(void)
{
    //debug
    int i=0;
    int count1;
    //
    Message message;
    MessageProcess(&gUART4ProtocolHandler);

    count1 = gUART4ProtocolHandler.txRawByteQueue.count;
    if (count1>0)
         flag_tx = 1;
#if 0
    if (RingBufferGetCount(&gUART4ProtocolHandler.txRawByteQueue))
    {
      //debug
      /*
        while  (RingBufferIsEmpty(&gUART4ProtocolHandler.txRawByteQueue) == false) 	//read for quere and pust to buf
        {
			uint8_t data;
			RingBufferDequeue(&gUART4ProtocolHandler.txRawByteQueue, &data);
			send_buf[i]=data;
            i++;

        }
        */
        //here:
        /*
        start senging flag

        */
        flag_tx = 1;

    }
#endif // 0

 #if 0  //edit by meng jun
    printf("count----->%d\n",count1);
    if(count1>0)
    {
            for (i=0;i<count1;i++)
                    printf(":%x-",send_buf[i]);
            printf("count----->end\n");
    }
 #endif

    if(MessageRead(&gUART4ProtocolHandler, &message)==ProtocolNoError)
    {
        #if 1
        //打印接收到的返回数据
        printf("Rx message: ");
        printf("message id:%d, rw:%d, isQueued:%d, paramsLen:%d\r\n",
                message.id, message.rw, message.isQueued, message.paramsLen);
        printf("params: ");
        for(int i=0; i<message.paramsLen; i++)
        {
            printf("%02x ", message.params[i]);
        }
        printf("\r\n");
        #endif
    }
}
