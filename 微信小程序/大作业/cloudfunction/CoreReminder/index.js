// 云函数入口文件
const cloud = require('wx-server-sdk');
cloud.init({ env: cloud.DYNAMIC_CURRENT_ENV }); // 使用当前云环境
const db = cloud.database();

// 云函数入口函数
exports.main = async (event, context) => {
  try {
    const res = await db.collection('reminder').where({
      done: false
    }).get();
    console.log('查询成功，返回的数据：', res.data);
    for (const item of res.data) {
      const timeflow = item.timeInput.split("-")
      const yearget=timeflow[0];
      const monthget=timeflow[1];
      const dayget=timeflow[2];
      const hourget=timeflow[3];
      const minuteget=timeflow[4];
      const date=new Date();
      const timeddl=`${yearget}-${monthget<10?`0${monthget}`:monthget}-${dayget<10?`0${dayget}`:dayget}T${hourget}:${minuteget}:00`;
      const timegetnow=`${date.getFullYear()}-${(date.getMonth()+1)<10?`0${date.getMonth()+1}`:date.getMonth()+1}-${date.getDate()<10?`0${date.getDate()}`:date.getDate()}T${((date.getHours()+8)%24)<10?`0${(date.getHours()+8)%24}`:(date.getHours()+8)%24}:${date.getMinutes()<10?`0${date.getMinutes()}`:date.getMinutes()}:00`;
      const dateset = new Date(timeddl);
      const datenow = new Date(timegetnow);
      const timeDifferenceMs = Math.abs(dateset - datenow);
      const timeDifferenceMinutes = Math.floor(timeDifferenceMs / (1000 * 60));
      const timestr=`${yearget}-${monthget<10?`0${monthget}`:monthget}-${dayget<10?`0${dayget}`:dayget} ${hourget}:${minuteget}:00`;
      let urgent_situation="紧急程度未知";
      switch(item.urgent_level) {
        case 1:
          urgent_situation = "不急";
          break;
        case 2:
          urgent_situation = "稍急";
          break;
        case 3:
          urgent_situation = "紧急";
          break;
        case 4:
          urgent_situation = "迫在眉睫";
          break;
        case 5:
          urgent_situation = "十万火急";
          break;
        default:
          urgent_situation = "无效的紧急程度等级";
          break;
      }
      if(timeDifferenceMinutes<=30){
      await db.collection('reminder').doc(item._id).update({
        data: {
          done: true,
        }
      });
      await cloud.callFunction({
        name: 'SendReminder',
        data:{
          titleValue:item.titleValue,
          urgent_level:urgent_situation,
          timeInput:timestr,
          _openid:item._openid
        }
      });
      
    }
    
  }
    return {
      success: true,
      message: '所有提醒已处理'
    };
  } catch (err) {
    console.error('处理出错：', err);
    return {
      success: false,
      message: '处理出错'
    };
  }
};