// 云函数入口文件
const cloud = require('wx-server-sdk');
cloud.init({ env: cloud.DYNAMIC_CURRENT_ENV }); // 使用当前云环境

// 云函数入口函数
exports.main = async (event, context) => {
  try {
        const result = await cloud.openapi.subscribeMessage.send({
          "touser":event._openid,
          "page": 'pages/detail_add/detail_add',
          "lang": 'zh_CN',
          "data": {
            "thing5": {
              "value": event.titleValue
            },
            "thing8": {
              "value":event.urgent_level
            },
            "thing2": {
              "value": "您的日程即将到期"
            },
            "date9": {
              "value": event.timeInput
            }
          },
          "templateId": 'v3J9x3kz6FYW3Xt4XwJ-d3qpbTU-J19XBh91HwOKqLA',
          "miniprogramState": 'developer',
        });
    
        return result;
      } catch (err) {
        console.error('发送订阅消息失败：', err);
        return {
          success: false,
          message: '发送订阅消息失败',
          error: err
        };
      }
}