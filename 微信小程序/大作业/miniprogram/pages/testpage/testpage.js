// pages/testpage/testpage.js

Page({

  /**
   * 页面的初始数据
   */
  data: {

  },

  /**
   * 生命周期函数--监听页面加载
   */
  onLoad(options) {
    const timeInput="2025-2-12-19-33"
    const timeflow = timeInput.split("-")
      const yearget=timeflow[0];
      const monthget=timeflow[1];
      const dayget=timeflow[2];
      const hourget=timeflow[3];
      const minuteget=timeflow[4];
      const date=new Date();
      const timeddl=`${yearget}-${monthget<10?`0${monthget}`:monthget}-${dayget<10?`0${dayget}`:dayget}T${hourget}:${minuteget}:00`;
      const timegetnow=`${date.getFullYear()}-${(date.getMonth()+1)<10?`0${date.getMonth()+1}`:date.getMonth()+1}-${date.getDate()<10?`0${date.getDate()}`:date.getDate()}T${date.getHours()<10?`0${date.getHours()}`:date.getHours()}:${date.getMinutes()<.10?`0${date.getMinutes()}`:date.getMinutes()}:00`;
    const dateset = new Date(timeddl);
      const datenow = new Date(timegetnow);
      const timeDifferenceMs = Math.abs(dateset - datenow);
      const timeDifferenceMinutes = Math.floor(timeDifferenceMs / (1000 * 60));
      console.log(timeddl);
      console.log(timegetnow);
    console.log(timeDifferenceMinutes);
    if(timeDifferenceMinutes<=30){
      console.log("wdf");
    }
  },

  /**
   * 生命周期函数--监听页面初次渲染完成
   */
  onReady() {

  },

  /**
   * 生命周期函数--监听页面显示
   */
  onShow() {

  },

  /**
   * 生命周期函数--监听页面隐藏
   */
  onHide() {

  },

  /**
   * 生命周期函数--监听页面卸载
   */
  onUnload() {

  },

  /**
   * 页面相关事件处理函数--监听用户下拉动作
   */
  onPullDownRefresh() {

  },

  /**
   * 页面上拉触底事件的处理函数
   */
  onReachBottom() {

  },

  /**
   * 用户点击右上角分享
   */
  onShareAppMessage() {

  }
})