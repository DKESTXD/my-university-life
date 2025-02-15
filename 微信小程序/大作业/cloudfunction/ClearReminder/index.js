// 云函数入口文件
const cloud = require('wx-server-sdk')

cloud.init({ env: cloud.DYNAMIC_CURRENT_ENV }) // 使用当前云环境
const db = cloud.database();
// 云函数入口函数
exports.main = async (event, context) => {
  try{
    await db.collection("reminder").where({
      done:true
    }).remove();

  }catch (err) {
    console.error('处理出错：', err);
    return {
      success: false,
      message: '处理出错'
    };
  }


}